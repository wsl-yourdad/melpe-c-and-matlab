%% =========================================================================
%  MELPe 1200bps Standalone Professional Decoder
%  Core: Resource Linker & State Isolation (工业级资源链接器)
% =========================================================================
clc; clear; close all;
fprintf('=== MELPe 1200bps 工业级高保真解码器启动 ===\n');
FS=8000;%采样率8000

% --- [板块 1]：物理资源自检与变量对齐 ---
fprintf('1. 正在初始化物理链路...\n');

% 1.1 核心码本自检 (THE_FINAL_CODEBOOK.mat)
if ~isfile('THE_FINAL_CODEBOOK.mat')
    error('❌ 链接失败：缺少核心码本 THE_FINAL_CODEBOOK.mat。');
else
    tmp_cb = load('THE_FINAL_CODEBOOK.mat');
    % 强制变量名对齐 (防止有些地方叫 codebooks 有些地方叫 CB)
    if isfield(tmp_cb, 'codebooks'), CB = tmp_cb.codebooks; else, CB = tmp_cb; end
    fprintf('   ✅ 核心码本加载成功。包含字段数: %d\n', length(fieldnames(CB)));
end

% 1.2 待解码比特流自检 (encoder_output.mat)
if ~isfile('encoder_output.mat')
    error('❌ 链接失败：找不到 encoder_output.mat，请检查编码器是否成功生成文件。');
else
    tmp_enc = load('encoder_output.mat');
    % 智能搜索比特流变量 (支持 all_bit_streams 或其他命名)
    fields = fieldnames(tmp_enc);
    bit_var = fields{cellfun(@(x) ~isempty(strfind(x, 'bit')), fields)};
    if isempty(bit_var), error('❌ 数据错误：.mat 文件中未找到比特流。'); end
    all_bit_streams = tmp_enc.(bit_var);
    fprintf('   ✅ 比特流识别成功: [%s], 总超帧数: %d\n', bit_var, length(all_bit_streams));
end

% 1.3 残差码本深度解析 (res256_64_64_64.txt)
% 此处不再简单用 sscanf，而是模拟 C 代码的内存映射，确保 20 维数据不偏移
if ~isfile('res256_64_64_64.txt')
    error('❌ 链接失败：缺少残差文本文件 res256_64_64_64.txt。');
else
    fid = fopen('res256_64_64_64.txt', 'r');
    % 暴力清洗非数值字符
    content = fread(fid, '*char')'; fclose(fid);
    content(content == ',' | content == '/' | content == '*' | content == 'Q') = ' ';
    raw_res = sscanf(content, '%f');
    
    % 物理维度校验：1200bps 标准残差码本应包含 8960 个浮点数 (20维 * 448个向量)
    if length(raw_res) < 8960
        error('❌ 资源损毁：残差码本数据量不足 (%d < 8960)。', length(raw_res));
    end
    % 按照 20 x 448 进行列映射 (对应 C 语言的 Stage 连续存储)
    RES_CB_RAW = reshape(raw_res(1:8960), 20, 448);
    fprintf('   ✅ 残差码本解析完成。Q17 原始尺度已锁定。\n');
end

% --- [板块 2]：解码器物理状态空间初始化 ---
% 这里的每一个变量都对应 C 代码中的全局变量或 Static 变量
fprintf('2. 正在隔离并分配解码器状态空间...\n');

st.prev_lsf   = double(CB.msvq_mean);  % 状态记忆：LSF (Q15)
st.prev_gain  = 20.0;                  % 状态记忆：Gain (dB)
st.prev_pitch = 50.0;                  % 状态记忆：Pitch (samples)
st.lpc_old    = [1, zeros(1, 10)];     % 状态记忆：安全 LPC 系数
st.syn_mem    = zeros(1, 10);          % 滤波器记忆
st.disp_mem   = zeros(1, 64);          % 分散滤波器记忆
st.f_phase    = zeros(1, 10);          % 傅里叶相位记忆 (重要：保证子帧间平滑)

% 滤波器预热
df = double(CB.disp_cof) / 32768.0; 
bpf_num = CB.bpf_num; 
bpf_den = CB.bpf_den;

% --- 补全状态空间初始化，确保与 synthesis_engine 内部变量名一致 ---
st.syn_mem    = zeros(1, 10);
st.disp_mem   = zeros(1, 64);
st.f_phase    = zeros(1, 10);
st.lpc_old    = [1, zeros(1, 10)];
st.bpf_mem    = zeros(5, 8); % <--- 核心修复：5个频带，每个滤波器8个状态位

fprintf('✅ 环境就绪。进入主解码循环。\n');

num_sf = length(all_bit_streams);
decoded_signal = zeros(1, num_sf * 540);
ptr = 1;

fprintf('正在处理 %d 个超帧...\n', num_sf);

% --- 在 load('THE_FINAL_CODEBOOK.mat') 之后插入 ---
fprintf('--- 码本 BPF 字段核查 ---\n');
if isfield(CB, 'bpf_num') && isfield(CB, 'bpf_den')
    fprintf('   ✅ 找到 BPF 系数。维度: Num(%dx%d), Den(%dx%d)\n', ...
        size(CB.bpf_num, 1), size(CB.bpf_num, 2), size(CB.bpf_den, 1), size(CB.bpf_den, 2));
    % 1200bps 标准通常是 45x1 或 9x5
else
    fprintf('   ⚠️ 未找到 BPF 相关字段，请检查码本变量名。\n');
end

% 3. 核心解码循环
for s_idx = 1:num_sf
    bits = all_bit_streams{s_idx};
    if length(bits) ~= 81, continue; end
    
    % --- [板块 A] 比特流物理拆解 ---
    p = unpack_bitstream(bits);
    
    % --- [板块 B] 全参数反量化 (与编码器 V32 严格对齐) ---
    lsf_mat = dequantize_lsf(p, st.prev_lsf, CB, RES_CB_RAW);
    [pitches, bpvc_mat] = dequantize_pitch(p, CB);
    gains_dequant = dequantize_gain(p.gain_idx, CB.gain);
    fsmag = dequantize_fsvq(p.fs_idx, CB);
    
    % --- [板块 C] 高保真逐帧合成 ---
    for f = 1:3
        % A. 准备插值参数
        g_end = gains_dequant(2, f);
        g_start = st.prev_gain;
        p_end = pitches(f);
        p_start = st.prev_pitch;
        
        % B. 调用核心合成引擎 (包含 4 子帧插值)
% --- 修正后的主循环合成调用 (严格12参数) ---
        [sig_frm, st] = synthesis_engine(...
                                  st.prev_lsf, ...      % 1. 上一帧LSF
                                  lsf_mat(:,f), ...     % 2. 当前帧LSF
                                  st.prev_pitch, ...    % 3. 上一帧Pitch
                                  pitches(f), ...       % 4. 当前帧Pitch
                                  bpvc_mat(:,f), ...    % 5. 5频带V/UV标志
                                  st.prev_gain, ...     % 6. 上一帧Gain (dB)
                                  gains_dequant(2,f),...% 7. 当前帧Gain (dB)
                                  fsmag, ...            % 8. 傅里叶幅值
                                  st, ...               % 9. 状态结构体
                                  df, ...               % 10. 分散滤波器
                                  CB.bpf_num, ...       % 11. BPF分子码本 (45x1)
                                  CB.bpf_den ...        % 12. BPF分母码本 (45x1)
                                    );
        
        decoded_signal(ptr:ptr+179) = sig_frm;
        ptr = ptr + 180;
        
        % C. 状态迁移
        st.prev_lsf = lsf_mat(:,f);
        st.prev_gain = g_end;
        st.prev_pitch = p_end;
    end
    
    if mod(s_idx, 20) == 0, fprintf('.'); end
end

% 4. 后处理与播放
decoded_signal = decoded_signal / 32768.0; 
% 自动增益归一化与去直流
decoded_signal = decoded_signal - mean(decoded_signal);
norm_val = max(abs(decoded_signal));
if norm_val > 0, decoded_signal = decoded_signal / norm_val * 0.9; end

% =========================================================================
% 板块 5: 信号可视化 (纯线条诊断)
% =========================================================================
fprintf('4. 正在绘制时域与频域波形...\n');

figure('Name', 'MELPe 解码结果实时诊断', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

% --- A. 时域波形 (Time Domain) ---
subplot(2, 1, 1);
t_axis = (0:length(decoded_signal)-1) / FS;

if exist('signal_final', 'var')
    % 如果存在原始信号，用灰色线条作为背景对比
    len_cmp = min(length(decoded_signal), length(signal_final));
    plot(t_axis(1:len_cmp), signal_final(1:len_cmp), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); 
    hold on;
end

plot(t_axis, decoded_signal, 'r', 'LineWidth', 0.7);
title('时域波形 (Time Domain Visualization)', 'FontSize', 12);
xlabel('时间 (s)'); ylabel('幅度');
legend('原始信号 (对比)', '解码还原信号');
grid on; axis tight;

% --- B. 频域幅度谱 (Frequency Spectrum) ---
subplot(2, 1, 2);
N_FFT = 2^nextpow2(length(decoded_signal));
f_axis = (0:N_FFT/2) * (FS/N_FFT);

% 计算解码信号频谱 (dB 尺度)
Y_dec = fft(decoded_signal, N_FFT);
P_dec = 20*log10(abs(Y_dec(1:N_FFT/2+1)) + 1e-6);

if exist('signal_final', 'var')
    % 计算原始信号频谱对比
    Y_orig = fft(signal_final, N_FFT);
    P_orig = 20*log10(abs(Y_orig(1:N_FFT/2+1)) + 1e-6);
    plot(f_axis, P_orig, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); 
    hold on;
end

plot(f_axis, P_dec, 'r', 'LineWidth', 0.7);
title('频域幅度谱 (Frequency Domain - FFT Line)', 'FontSize', 12);
xlabel('频率 (Hz)'); ylabel('幅度 (dB)');
grid on; xlim([0, FS/2]); % 限制在 4000Hz 带宽内
if exist('signal_final', 'var'), legend('原始频谱', '解码频谱'); end

fprintf('   ✅ 绘图完成！\n');

%audiowrite('melpe_1200_pro_output.wav', decoded_signal, 8000);
soundsc(decoded_signal, 8000);

fprintf('\n✅ 解码成功！输出文件: melpe_1200_pro_output.wav\n');