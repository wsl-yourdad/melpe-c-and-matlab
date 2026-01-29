%% =========================================================================
%  MELPe 1200bps 语音编码仿真 - 模块 I: 初始化与加载 (优化版)
%  
%  功能: 
%    1. 加载码本与配置参数
%    2. 读取语音并进行预处理 (60Hz高通)
%    3. 严格按照超帧 (3帧/组) 结构进行对齐裁剪
% =========================================================================

clc; clear; close all;
fprintf('=== MELPe 1200bps 仿真启动 ===\n');

%% 一。加载码本数据库
% -------------------------------------------------------------------------
codebook_file = 'THE_FINAL_CODEBOOK.mat';

if ~isfile(codebook_file)
    error('❌ 错误: 找不到 %s。请确保你已经运行了数据构建脚本！', codebook_file);
end

fprintf('正在加载码本...');
load(codebook_file); % 加载得到结构体 'codebooks'

% === V31 修复: 正确加载 20维 Q17 残差码本
% ===%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- 加载基础码本 ---
load('THE_FINAL_CODEBOOK.mat'); % 此时产生了 codebooks 变量


fprintf(' 完成。\n');

% --- 完整性自检 ---
%required_fields = {'inpCoef', 'pitch_vvv', 'gain', 'lsp_vq', 'res_vq', 'fec_dec'};
%missing = setdiff(required_fields, fieldnames(codebooks));
%if ~isempty(missing)
%    error('❌ 码本数据不完整！缺少字段: %s', strjoin(missing, ', '));
%end
%fprintf('✅ 码本完整性校验通过。\n');

% 提取滤波器系数
LPF_NUM = codebooks.lpf_num; LPF_DEN = codebooks.lpf_den;
BPF_NUM = codebooks.bpf_num; BPF_DEN = codebooks.bpf_den;

%% 打印并校验码本内容
% -------------------------------------------------------------------------
fprintf('\n--- 码本详细信息查询 ---\n');
fields = fieldnames(codebooks);
fprintf('包含字段数: %d\n', length(fields));

% 创建表格显示字段名和数据维度，方便查阅
fprintf('%-15s | %-15s | %-10s\n', '字段名', '维度', '数据类型');
fprintf('%s\n', repmat('-', 1, 45));

for i = 1:length(fields)
    val = codebooks.(fields{i});
    sz = size(val);
    sz_str = sprintf('%dx%d', sz(1), sz(2));
    fprintf('%-15s | %-15s | %-10s\n', fields{i}, sz_str, class(val));
end

% --- 核心量化表逻辑抽检 (增强版) ---
if isfield(codebooks, 'pitch_vvv')
    fprintf('\n[基音码本抽检] (对应附录表索引):\n');
    
    % 确保索引不越界
    check_idx = [1, 17, 99]; 
    max_len = length(codebooks.pitch_vvv);
    
    for idx = check_idx
        if idx <= max_len
            val = codebooks.pitch_vvv(idx);
            
            % 检查数值是否有效（防止 NaN 或 Inf 导致 dec2hex 报错）
            if isfinite(val)
                % 强制转换为 uint32 以确保 dec2hex 能够处理
                hex_str = dec2hex(uint32(val));
                fprintf('  索引 %2d -> 码字: 0x%s\n', idx-1, hex_str);
            else
                fprintf('  索引 %2d -> 码字: 数据无效(NaN/Inf)\n', idx-1);
            end
        end
    end
end

% 检查 1200bps 特有的 LSP 矢量量化维度
if isfield(codebooks, 'lsp_vq')
    % 1200bps 40bits 通常由多级 VQ 组成，检查其行数
    fprintf('[LSP 码本检测]: 维度为 %dx%d\n', size(codebooks.lsp_vq,1), size(codebooks.lsp_vq,2));
end
fprintf('%s\n', repmat('-', 1, 45));

%% 二. 核心参数配置
% -------------------------------------------------------------------------
FS = 8000;               % 采样率 8kHz
FRAME_LEN = 180;         % 单帧长度 (22.5ms)
SUPERFRAME_SIZE = 3;     % 超帧包含 3 帧 (67.5ms)
LPC_ORD = 10;            % LPC 阶数

fprintf('\n--- 参数配置 ---\n');
fprintf('采样率 (FS):       %d Hz\n', FS);
fprintf('超帧结构:          %d 帧/组 (%.1f ms)\n', SUPERFRAME_SIZE, SUPERFRAME_SIZE*FRAME_LEN/FS*1000);


%% 三. 读取语音文件
% -------------------------------------------------------------------------
wav_file = 'input_8k_5s.wav'; % 建议使用你的文件名

try
    [sig_raw, actual_fs] = audioread(wav_file);
    
    % 采样率转换
    if actual_fs ~= FS
        sig_raw = resample(sig_raw, FS, actual_fs);
    end
    
    % 单声道化
    if size(sig_raw, 2) > 1
        sig_raw = sig_raw(:, 1);
    end
    
    % 截取前5秒进行仿真
    points_5s = FS * 5;
    signal_5s = sig_raw(1:min(points_5s, length(sig_raw)));
    fprintf('✅ 音频加载成功，截取前5秒分析。\n');
catch
    error('❌ 读取音频失败，请检查文件名！');
end


%% 四。 预处理与超帧对齐 (改进部分)
% -------------------------------------------------------------------------
% 4.1 执行预处理 (60Hz高通滤波，去除直流和低频噪声)
% 建议在裁剪前滤波，以消除滤波器暂态响应对边界的影响
fprintf('执行预处理 (高通滤波)... ');
filtered_signal_full = melp_preprocessing(signal_5s, FS);
fprintf('完成。\n');

% 4.2 执行加扰实验 (仅用于论文复现对比，不参与后续编码)
disturbed_signal = add_interference(signal_5s, FS);

% 4.3 超帧对齐裁剪
% 计算总样点数对应的最大整数超帧数
total_samples = length(filtered_signal_full);
samples_per_superframe = FRAME_LEN * SUPERFRAME_SIZE; % 180 * 3 = 540 点 一个超帧对应的点数

num_superframes = floor(total_samples / samples_per_superframe); %5s样本内超帧个数
num_frames = num_superframes * SUPERFRAME_SIZE;%5s样本里的总帧数

% 对信号进行最终裁剪，使其严格符合超帧边界
final_sample_count = num_superframes * samples_per_superframe; %超帧个数*540（一个超帧对应的点数）
signal_final = filtered_signal_full(1:final_sample_count); % 后续模块使用此变量

fprintf('--- 超帧对齐完成 ---\n');
fprintf('  -> 总超帧数: %d\n', num_superframes);
fprintf('  -> 总帧数:   %d\n', num_frames);
fprintf('  -> 最终时长: %.2f 秒\n', length(signal_final)/FS);


% 4.4 可视化结果 (复现论文图 2.5)
% -------------------------------------------------------------------------
figure('Name', 'MELPe 预处理与对齐监控', 'Position', [100 50 900 700]);

subplot(3,1,1);
plot(signal_5s); title('① 原始语音信号 (5秒截取)');
grid on; ylabel('幅度');

subplot(3,1,2);
plot(disturbed_signal); title('② 加扰实验信号 (50Hz+直流)');
grid on; ylabel('幅度');

subplot(3,1,3);
plot(signal_final); title('③ 预处理并超帧对齐后的信号 (用于编码器输入)');
grid on; ylabel('幅度'); xlabel('采样点');

fprintf('✅ 模块 I 初始化完成。当前信号已准备好进入基音提取和 LPC 分析。\n');
%%五。 加入的调试步骤
fprintf('--- 信号能量检查 ---\n');
fprintf('  signal_final 最大幅值: %.4f\n', max(abs(signal_final)));
fprintf('  signal_final 平均能量: %.4f\n', mean(signal_final.^2));
if max(abs(signal_final)) < 0.01
    warning('⚠️ 信号能量太弱，可能导致基音提取失败！');
    % 尝试进行归一化
    signal_final = signal_final / max(abs(signal_final));
    fprintf('  ✅ 已自动执行归一化处理。\n');
end


%% =========================================================================
%  六。 核心编码循环 (已修复逻辑错误的完美版)
% =========================================================================
% --- 初始化编码器状态与C标准容器 ---
fprintf('--- 初始化编码器状态与C标准容器 ---\n');

% [关键修复] 初始化 LSF 预测器记忆
if isfield(codebooks, 'msvq_mean')
    last_lsf_quant = codebooks.msvq_mean; 
else
    error('❌ 码本中缺少 msvq_mean，无法初始化 LSF 预测器！');
end

num_total_frames = num_superframes * 3;
% 初始化可视化容器
viz_pitches = zeros(1, num_total_frames);
viz_voicing = zeros(1, num_total_frames);
viz_gains   = zeros(1, num_total_frames);
viz_fsmags  = zeros(10, num_total_frames);
viz_p1      = zeros(1, num_total_frames); % 调试用
viz_p2      = zeros(1, num_total_frames); % 调试用
viz_g1      = zeros(1, num_total_frames); % 调试用

% 初始化比特流容器
superframe_lsf_indices_cell = cell(1, num_superframes);
superframe_pitch_idx_c      = zeros(1, num_superframes);
superframe_bpvc_idx_c       = zeros(1, num_superframes);
superframe_jitter_idx_c     = zeros(1, num_superframes);
superframe_gain_indices     = zeros(1, num_superframes);
all_bit_streams             = cell(1, num_superframes);

% ====================[ 主循环开始 ]====================
for s_idx = 1:num_superframes
    
    % --- Step 1: 子帧特征提取 (循环收集3帧数据) ---
    sf_lsfs         = zeros(LPC_ORD, 3);
    sf_pitches      = zeros(1, 3);
    sf_v_dec        = zeros(1, 3);
    sf_gains        = zeros(2, 3);
    sf_fsmags_tmp   = zeros(10, 3);
    lpc_a_superframe = zeros(LPC_ORD, 3);
    sf_bpvc         = zeros(5, 3);
% --- Step 1: 子帧特征提取 (最终修正版：修复维度+数值爆炸) ---
    for f_in_s = 1:3
        frame_idx = (s_idx-1) * 3 + f_in_s;
        start_pt = (frame_idx-1) * FRAME_LEN + 1;
        curr_frame = signal_final(start_pt : start_pt + FRAME_LEN - 1);
        
        % =================================================================
        % 1. LPC 分析 & 残差计算
        % =================================================================
        win = hamming(FRAME_LEN);
        sp_win = curr_frame .* win;
        
        % 计算自相关
        r_full = xcorr(sp_win, LPC_ORD, 'biased');
        r = r_full(LPC_ORD+1:end); 
        r(1) = r(1) * 1.0001; % 防止奇异
        
        % Levinson 求解
        a = levinson(r, LPC_ORD); 
        
        % [BWEX] 带宽扩展
        GAMMA = 0.994;
        bw_exp = GAMMA .^ (0:LPC_ORD);
        a_bw = a .* bw_exp; 
        
        % [存入矩阵] 修复维度问题：取后10个，并转置为列向量
        lpc_a_superframe(:, f_in_s) = a_bw(2:end).'; 
        
        % 计算残差 (用于增益计算)
        res_frame = filter(a_bw, 1, curr_frame);
        
        % =================================================================
        % 2. LSF 计算 (使用 MATLAB 自带函数 poly2lsf 确保是弧度!)
        % =================================================================
        % poly2lsf 返回的是弧度 0 ~ pi
        lsf_rad = poly2lsf(a_bw); 
        
        % 映射到 Q15 (0~pi -> 0~32768)
        % 只有当 lsf_rad 是弧度时，这个公式才对！
        lsf_q15_int = (lsf_rad / pi) * 32768.0;
        
        % 钳位 (0~32767)
        lsf_q15_int = sort(lsf_q15_int);
        lsf_q15_int(lsf_q15_int < 50) = 50;
        lsf_q15_int(lsf_q15_int > 32700) = 32700;
        
        % [修复报错] 强制转为列向量赋值，防止维度不兼容
        sf_lsfs(:, f_in_s) = lsf_q15_int(:);
        
        % =================================================================
        % 3. 基音与增益 (保持不变)
        % =================================================================
        [P1, ~] = melp_pitch_integer(curr_frame, FS);
        [P2, Vp2] = melp_pitch_fraction(curr_frame, P1, FS);
        [P3, ~, Voicing_flag] = melp_pitch_final(curr_frame, FS, P2, Vp2);
        
        if Voicing_flag == 1 && P3 < 10, P3 = max(P2, 20); end
        
        sf_pitches(f_in_s) = P3; 
        sf_v_dec(f_in_s) = Voicing_flag;
        
        % 增益输入放大以匹配C代码量级
        [G1, G2] = melp_gain_calculator(curr_frame*32767, res_frame*32767, P3, Voicing_flag);
        sf_gains(:, f_in_s) = [G1; G2];
        
        % 4. 频带分析 (修改版：计算真实 BPVC)
        % =================================================================
        if Voicing_flag == 1
            [fsmag_norm, ~] = melp_harmonic_magnitudes(a_bw, P3, FS);
            sf_fsmags_tmp(:, f_in_s) = fsmag_norm;
            
            % [新增] 计算当前帧的 5-band BPVC
            % 必须传入 codebooks 以获取滤波器系数
            curr_bpvc = melp_bpvc_analysis(curr_frame, FS, codebooks,P3);
        else
            sf_fsmags_tmp(:, f_in_s) = ones(10, 1);
            
            % [新增] 清音帧 BPVC 全为 0
            curr_bpvc = zeros(5, 1);
        end
        
        % [存入容器]
        sf_bpvc(:, f_in_s) = curr_bpvc;
        
        % 调试数据
        viz_pitches(frame_idx) = P3;
        viz_voicing(frame_idx) = Voicing_flag;
        viz_gains(frame_idx)   = G2;
    end
    
    % =====================================================================
    % Step 2: 超帧联合量化 (只有当3帧数据都齐了才执行！)
    % =====================================================================
    
    % --- 2.1 Pitch, BPVC, Jitter ---
    lpc_coeffs_full = [ones(1, 3); lpc_a_superframe]; % 凑成 11x3 矩阵
% 原来的调用：
% [p_idx_c, bpvc_idx_c, jit_idx_c, uv_bits_c] = quantize_pitch_bpvc_c_style(sf_pitches, sf_v_dec, lpc_coeffs_full, codebooks);

% ✅ 新的调用
% **************************************************************************************************************************************************************************
    % [修正] 传入第 3 帧的真实 BPVC 数据，不再传空值 []
    bpvc_target = sf_bpvc(:, 3); 
    [p_idx_c, bpvc_idx_c, jit_idx_c,encoded_uv_c] = quantize_pitch_bpvc_c_style(sf_pitches, sf_v_dec, bpvc_target,codebooks);

% 提示：uv_bits_c 不需要了，因为 packer 会自己生成。
        
    superframe_pitch_idx_c(s_idx)  = p_idx_c;
    superframe_bpvc_idx_c(s_idx)   = bpvc_idx_c;
    superframe_jitter_idx_c(s_idx) = jit_idx_c;
    
    % --- 2.2 Gain (10-bit VQ) ---
    target_gain_vec_6D = reshape(sf_gains, 6, 1);
    superframe_gain_indices(s_idx) = gain_vq_search(target_gain_vec_6D, codebooks.gain);
    
    % --- 2.3 LSF 量化 (最关键的修复) ---
    % 注意：必须接收第二个返回值 lsf_quant_out 用于更新记忆
    % === 🔍 间谍代码：查看数值范围 ===
     fprintf('\n[DEBUG] LSF 数值范围诊断:\n');
     fprintf('  -> 输入 LSF (Hz) 前5个值: %.2f %.2f %.2f %.2f %.2f\n', sf_lsfs(1:5, 3));
     fprintf('  -> 码本 (9bit) 前5个值: %.4f %.4f %.4f %.4f %.4f\n', codebooks.lsp_uv_9(1:5, 1));
     fprintf('  -> 码本最大值: %.4f\n', max(codebooks.lsp_uv_9(:)));
% ================================
    [lsf_indices_c, lsf_quant_out] = quantize_lsf_c_style(sf_lsfs, sf_v_dec, codebooks, last_lsf_quant);
    superframe_lsf_indices_cell{s_idx} = lsf_indices_c;
    
    % [重要] 更新 LSF 预测器的状态记忆！没有这一步 LSF 就死锁了！
    last_lsf_quant = lsf_quant_out; 

    %检测LSF部分
% === LSF 物理合理性验证 (修复版) ===
fprintf('\n[物理检查] 超帧 %d 的 LSF 状态更新值:\n', s_idx);

% 检查 lsf_quant_out 的维度
[rows, cols] = size(lsf_quant_out);

    if cols == 3
    % 如果函数返回了全部三帧
        for f = 1:3
        val = lsf_quant_out(:, f);
        check_lsf_values(val, f);
        end
    else
    % 如果函数只返回了 F3 (当前 V32.1 的逻辑)
    val = lsf_quant_out; 
    check_lsf_values(val, 3); % 标记为第3帧
end
    
    % --- 2.4 Fourier Series ---
    fourier_vq_idx = 1;
    last_voiced_frame_idx = find(sf_v_dec == 1, 1, 'last');
    if ~isempty(last_voiced_frame_idx)
         fsmag_target_vec = sf_fsmags_tmp(:, last_voiced_frame_idx);
         if isfield(codebooks, 'fsvq_cb')
            fourier_vq_idx = melp_vq_pitch_search(fsmag_target_vec, codebooks.fsvq_cb);
         end
    end
    
    % --- Step 3: 最终打包 ---
% ✅ 修正后的正确调用：
% =====================================================================
    % Step 3: 最终打包 (V33 修正调用顺序)
    % 必须严格对应 pack_bits_c_style 的输入参数顺序:
    % (pitch, gain, uv, lsf_indices, jitter, bpvc, fs)
    % =====================================================================
    
    % 准备参数
    p_idx = p_idx_c; % Pitch
    g_idx = superframe_gain_indices(s_idx); % Gain%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lsf_idx_struct = superframe_lsf_indices_cell{s_idx}; % LSF Struct
    jit_idx = superframe_jitter_idx_c(s_idx); % Jitter
    bpvc_idx = superframe_bpvc_idx_c(s_idx); % BPVC
    fs_idx = fourier_vq_idx - 1; % FSVQ (注意 0-based)
    
   % =====================================================================*****************************************************
    % [调试埋点] 编码器真值打印 (修正版)
    % =====================================================================
    % 智能获取 LSF 第一个数值，不管字段名是 f1_s1 还是 idx1
    try
        % 1. 获取结构体字段名列表
        lsf_fields = fieldnames(lsf_idx_struct);
        % 2. 取第一个字段的名字
        first_field = lsf_fields{1}; 
        % 3. 取第1帧的该字段数值
        if length(lsf_idx_struct) > 1
            lsf_val = lsf_idx_struct(1).(first_field);
        else
            lsf_val = lsf_idx_struct.(first_field);
        end
        lsf_debug_tag = sprintf('LSF_1:%-3d', lsf_val);
    catch
        lsf_debug_tag = 'LSF:ERR';
    end
    
    % 打印核心参数 (格式: SF号 | UV | P | G | LSF | BPVC | J)
    fprintf('[ENC_TRUTH] SF:%03d | UV:%s | P:%03d | G:%04d | %-10s | BPVC:%02d | J:%d\n', ...
        s_idx, ...
        sprintf('%d%d%d', encoded_uv_c), ...
        p_idx, ...
        g_idx, ...
        lsf_debug_tag, ...
        bpvc_idx, ...
        jit_idx);
    % =====================================================================*******************************************

    % 调用打包函数
    all_bit_streams{s_idx} = pack_bits_c_style(...
        p_idx, ...      % 1. Pitch
        g_idx, ...      % 2. Gain
        encoded_uv_c, ...  % 3. UV
        lsf_idx_struct, ... % 4. LSF Indices (Struct)
        jit_idx, ...    % 5. Jitter
        bpvc_idx, ...   % 6. BPVC
        fs_idx ...      % 7. FSVQ
    );
    
    % 简单的进度条
    if mod(s_idx, 10) == 0
        fprintf('已处理超帧: %d / %d\n', s_idx, num_superframes);
    end
    
end
fprintf('✅ 核心编码循环结束！\n');
%% =========================================================================
%  6. 仿真结果可视化 (超帧编码特性分析)
%  功能: 对比原始信号与提取的特征参数 (基音、增益、清浊音)
% =========================================================================
fprintf('\n正在生成 1200bps 编码特征可视化图表...\n');

figure('Name', 'MELPe 1200bps 编码特征监控', 'Position', [100 50 1100 700]);% 左下宽高

% --- 1. 时域信号对比 (带超帧边界线) ---
subplot(4,1,1);
t_axis = (0:length(signal_final)-1) / FS;
plot(t_axis, signal_final, 'Color', [0.5 0.5 0.5]); hold on;
% 绘制超帧边界 (每 67.5ms 一道红虚线)
sf_boundary = (0:num_superframes) * (SUPERFRAME_SIZE * FRAME_LEN / FS);%计算超帧边界 i*3*180/8000
for b = 1:length(sf_boundary)
    % line([x1 x2], [y1 y2])：在点(x1, y1)和(x2, y2)之间画线。
    line([sf_boundary(b) sf_boundary(b)], [-1 1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
end
title('① 原始语音信号 (红色虚线为 67.5ms 超帧边界)');
ylabel('幅度'); grid on; xlim([0 t_axis(end)]);

% --- 2. 基音周期轨迹 (跨子帧) ---
subplot(4,1,2);
frame_t = (0:num_total_frames-1) * (FRAME_LEN / FS) + (FRAME_LEN/FS/2);
% 仅绘制浊音部分的基音
pitch_for_plot = viz_pitches;
pitch_for_plot(viz_voicing == 0) = NaN; 
plot(frame_t, pitch_for_plot, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
title('② 基音周期轨迹 (仅显示浊音帧)');
ylabel('周期 (采样点)'); grid on; xlim([0 t_axis(end)]); ylim([0 180]);

% % --- 3. 增益包络 (G2 语音增益) ---
subplot(4,1,3);
plot(frame_t, viz_gains, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 4);
title('③ 增益包络演变 (G2 Gain)');
ylabel('分贝 (dB)'); grid on; 
xlim([0 t_axis(end)]);
% 允许 Y 轴根据实际分贝值(通常在 10~80dB)自适应
axis tight;

% --- 3. 清浊音决策 (V/UV) ---
subplot(4,1,4);
stem(frame_t, viz_voicing, 'm', 'LineWidth', 1.2);
title('④ 清浊音判决标志 (1=浊音, 0=清音)');
ylabel('Voicing'); xlabel('时间 (s)'); grid on;
xlim([0 t_axis(end)]); ylim([-0.2 1.2]);

% 整体大标题
sgtitle(['MELPe 1200bps 仿真分析 - 共 ', num2str(num_superframes), ' 组超帧'], ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('✅ 可视化完成！请检查图形窗口以分析编码特征。\n');


%% =========================================================================
%% =========================================================================
%  7. [新] 命令行数据抽检 (C标准兼容版)
% =========================================================================
fprintf('\n=== [C-Standard 核心数据抽检] 第 1-20 帧特征提取明细 ===\n');
fprintf('帧号\tP3(终)\t清浊\t增益G2(dB)\tLSF索引(前2)\tPitch索引\tBPVC索引\tJitter\n');
fprintf('--------------------------------------------------------------------------------------------------\n');

start_f = 1;
end_f = min(20, num_total_frames);

for idx = start_f : end_f
    s_idx = floor((idx-1)/3) + 1; % 计算当前帧所属的超帧
    
    % 安全地从cell中提取LSF索引
    lsf_indices_str = 'N/A';
% --- melp_1200bps_2.m 主循环末尾 ---
% 替换掉原来的那行 sprintf
    if ~isempty(lsf_indices_c(1).idx1)
    lsf_info = sprintf('M1_F1Idx:%d', lsf_indices_c(1).idx1);
    else
    lsf_info = sprintf('M2_R1Idx:%d', lsf_indices_c(1).res_idx);
    end
% 然后再进行 fprintf 打印

    fprintf('%d\t%.2f\t%d\t[%.1f]\t\t%-12s\t%-9d\t%-8d\t%d\n', ...
            idx, ...
            viz_pitches(idx), ...
            viz_voicing(idx), ...
            viz_gains(idx), ...
            lsf_indices_str, ...
            superframe_pitch_idx_c(s_idx), ...
            superframe_bpvc_idx_c(s_idx), ...
            superframe_jitter_idx_c(s_idx));
end
fprintf('--------------------------------------------------------------------------------------------------\n');

%% =========================================================================
% 8. [新] 81-bit 码流字段拆解抽检 (C标准)
% =========================================================================
fprintf('\n=== [81-bit C-Standard 码流字段拆解抽检] ===\n');
fprintf('%-6s | %-1s | %-4s | %-9s | %-42s | %-10s | %-6s | %-8s | %-1s\n', ...
'超帧', 'S', 'UV+P', 'Pitch(9)', 'LSF(42)', 'Gain(10)', 'BPVC(6)', 'FS(8)', 'J');
fprintf('%s\n', repmat('-', 1, 95));

for s_idx = 1:num_superframes
    if mod(s_idx, 10) == 0 || s_idx == 1
        if s_idx <= length(all_bit_streams) && ~isempty(all_bit_streams{s_idx})
            b_str = all_bit_streams{s_idx};
            
            % 严格按照 pack_bits_c_style 的顺序进行切分
            f_sync   = b_str(1);
            f_uv     = b_str(2:5);
            f_pitch  = b_str(6:14);
            f_lsf    = b_str(15:56);
            f_gain   = b_str(57:66);
            f_bpvc   = b_str(67:72);
            f_fs     = b_str(73:80);
            f_jitter = b_str(81);
            
            fprintf('SF %2d  | %s | %s | %s | %s... | %s | %s | %s | %s\n', ...
                    s_idx, f_sync, f_uv, f_pitch, f_lsf(1:15), f_gain, f_bpvc, f_fs, f_jitter);
        end
    end
end
fprintf('%s\n', repmat('-', 1, 95));

%% =========================================================================
%  10. [核心验证] LSF "灵魂" 回环测试 (Residue Resynthesis)
%  原理: 使用 量化后的LSF + 原始残差。
%  预期: 如果LSF量化正确，声音应该极其清晰，几乎无损，只有轻微的量化底噪。
%        如果还有爆鸣/蟋蟀音，说明 LSF 量化/反量化 依然有错。
% =========================================================================
fprintf('\n=== 正在执行 LSF 灵魂验证 (Residual Resynthesis) ===\n');

% 准备合成容器
synth_debug = zeros(size(signal_final));
lpc_last_debug = zeros(10,1); 

% 1. 全局计算原始残差 (Ideal Excitation)
% 为了验证方便，我们重新用未量化的 LPC 算一遍最纯净的残差
% 也可以直接用 signal_final，在循环里 filter
resid_full = zeros(size(signal_final));

fprintf('正在合成...\n');

for s_idx = 1:num_superframes
    
    % --- A. 解码 LSF (模拟解码器行为) ---
    % 从编码器输出的索引中恢复 LSF
    indices = superframe_lsf_indices_cell{s_idx};
    uv_flag = [str2double(all_bit_streams{s_idx}(2)), ...
               str2double(all_bit_streams{s_idx}(3)), ...
               str2double(all_bit_streams{s_idx}(4))];
               
    % 这里我们需要一个简易的解码函数 (或者手动查表)
    % 为简单起见，直接调用你写好的 quantize 函数的"反量化"逻辑
    % 但更真实的是直接用本次循环产生的 sf_lsfs (如果是量化后的)
    
    % 咱们直接用本次编码循环里计算出的 "sf_lsfs" (这是量化前的) 
    % 或者如果有量化后的值更好。
    % 这里为了测试"量化器"是否坏了，我们需要用"量化后"的值。
    % 但你的主循环里没存量化后的 LSF 向量，只存了索引。
    
    % 没关系，我们先验证 LSF 转换逻辑本身对不对！
    % 我们直接用 sf_lsfs (Q15) -> 还原为 LPC -> 合成
    
    % 为了严谨，我们这里重新进行简单的 LPC 分析 -> LSF -> LPC 转换
    % 看看这一套流程有没有引入失真
    
    for f = 1:3
        idx_start = (s_idx-1)*3*FRAME_LEN + (f-1)*FRAME_LEN + 1;
        curr_frm = signal_final(idx_start : idx_start + FRAME_LEN - 1);
        
        % 1. 重新分析 (Analysis)
        win = hamming(FRAME_LEN);
        r = xcorr(curr_frm .* win, LPC_ORD, 'biased');
        a_ana = levinson(r(LPC_ORD+1:end), LPC_ORD);
        
        % [关键] 加带宽扩展 (跟编码器一致)
        bw = 0.994 .^ (0:10);
        a_bw = a_ana .* bw;
        
        % 2. 计算残差 (Residual) = 原始语音 * A(z)
        % 这就是完美的激励源！
        res_frm = filter(a_bw, 1, curr_frm);
        
        % 3. LSF 转换回环测试 (LPC -> LSF -> LPC)
        lsf_rad = poly2lsf(a_bw);
        
        % [模拟量化过程] 转Q15再转回来
        lsf_q15 = (lsf_rad / pi) * 32768.0;
        lsf_q15 = round(lsf_q15); % 模拟整数截断
        
        % [反向还原] Q15 -> 弧度 -> LPC
        lsf_rec_rad = (lsf_q15 / 32768.0) * pi;
        a_rec = lsf2poly(lsf_rec_rad);
        
        % 4. 合成 (Synthesis) = 残差 / A_rec(z)
        [syn_frm, lpc_last_debug] = filter(1, a_rec, res_frm, lpc_last_debug);
        
        synth_debug(idx_start : idx_start + FRAME_LEN - 1) = syn_frm;
    end
end

% 归一化并播放
synth_debug = synth_debug / max(abs(synth_debug)) * 0.95;
audiowrite('debug_residual_resynth.wav', synth_debug, FS);

soundsc(double(real('debug_residual_resynth.wav')), FS);

fprintf('✅ 验证文件已生成: debug_residual_resynth.wav\n');
fprintf('   -> 请听这个文件！\n');
fprintf('   -> 如果清晰: 说明 LPC/LSF 转换、BWEX、Q15映射 全都对了！\n');
fprintf('   -> 如果有爆鸣: 说明 poly2lsf/lsf2poly 或者 Q15 映射还有问题。\n');

% 画图对比
figure('Name', 'LSF Soul Check');
subplot(2,1,1); plot(signal_final(1:1000)); title('Original'); grid on;
subplot(2,1,2); plot(synth_debug(1:1000)); title('Resynthesized (LPC Loopback)'); grid on;
% =========================================================================
% 9. [最终] 保存编码器输出以供解码器使用
% =========================================================================
fprintf('\n正在为解码器保存输出...\n');
save('encoder_output_new.mat', 'all_bit_streams', 'signal_final', 'FS');
fprintf('✅ encoder_output_new.mat 已保存。\n');
% =========================================================================
% [新增] LSF 钳位函数 (C-Style lpc_clmp)
% 作用：强制 LSF 排序并保持最小间隔，防止解码器滤波器爆炸
% =========================================================================
function lsf = lpc_clamp(lsf, delta)
    ORDER = length(lsf);
    
    % 1. 强制排序 (Sort)
    lsf = sort(lsf);
    
    % 2. 强制最小间隔 (Ensure minimum separation)
    % C代码逻辑复刻：循环10次确保间隔
    for iter = 1:10
        unsorted = false;
        for i = 1:ORDER-1
            if lsf(i+1) - lsf(i) < delta
                mid = (lsf(i+1) + lsf(i)) / 2;
                lsf(i)   = mid - delta/2;
                lsf(i+1) = mid + delta/2;
                unsorted = true;
            end
        end
        % 边界检查 (0 ~ 0.5 Normalized Hz)
        if lsf(1) < delta, lsf(1) = delta; end
        if lsf(end) > 0.5-delta, lsf(end) = 0.5-delta; end
        
        if ~unsorted, break; end
    end
end

function [best_idx, q_vec] = vq_search_general(target, codebook)
    % 无论输入是 10 维还是 20 维，都能正常计算欧氏距离
    num_entries = size(codebook, 2);
    min_dist = inf;
    best_idx = 1;
    for i = 1:num_entries
        diff = target - codebook(:, i);
        dist = sum(diff.^2); % 距离计算
        if dist < min_dist
            min_dist = dist;
            best_idx = i;
        end
    end
    q_vec = codebook(:, best_idx);
end

% 辅助本地函数 (可以直接放在循环里或脚本末尾)
function check_lsf_values(val, frame_idx)
    fprintf('  帧 %d: ', frame_idx);
    if any(val < 0 | val > 32768)
        fprintf('❌ 溢出! 值范围: [%.1f, %.1f]\n', min(val), max(val));
    elseif any(diff(val) < 0)
        fprintf('⚠️ 非单调! 滤波器必炸。\n');
    else
        fprintf('✅ 正常 (范围: %.1f - %.1f)\n', min(val), max(val));
    end
end

function bpvc = melp_bpvc_analysis(sig, FS, cb, pitch)
% MELP 5频带清浊音分析器 (High Fidelity Version)
% 输入: 
%   sig: 当前语音帧 (通常 180 点)
%   fs: 采样率 (8000)
%   cb: 码本结构体 (含滤波器系数)
%   pitch: 当前帧的基音周期 (P3) - 【新增关键参数】
% 输出: 
%   bpvc: 5x1 向量 (1=Voiced, 0=Unvoiced)

    bpvc = zeros(5, 1);
    
    % 1. 基础检查
    if ~isfield(cb, 'bpf_num') || ~isfield(cb, 'bpf_den')
        bpvc = ones(5, 1); return; % 容错
    end
    
    % 2. 准备数据
    % 解析滤波器 (45x1 -> 9x5)
    b_num = reshape(cb.bpf_num, 9, 5);
    b_den = reshape(cb.bpf_den, 9, 5);
    
    % 确定相关性计算的滞后点数 (Lag)
    % pitch 是浮点数，转为整数 lag
    pitch_lag = round(pitch);
    if pitch_lag < 20, pitch_lag = 20; end
    if pitch_lag > 160, pitch_lag = 160; end
    
    % 3. 逐频带分析
    for k = 1:5
        % A. 带通滤波
        % 注意：为了消除滤波器暂态，理想情况下应维护滤波器状态 (zi)
        % 这里为了简化且不破坏主循环结构，直接滤波 (对于短帧稍有误差但可接受)
        band_sig = filter(b_num(:, k), b_den(:, k), sig);
        
        % B. 计算 Pitch-Lag 归一化自相关 (核心!)
        % 公式: R = sum(x[n]*x[n-T]) / sqrt(sum(x[n]^2)*sum(x[n-T]^2))
        
        N = length(band_sig);
        start_idx = pitch_lag + 1;
        
        if start_idx >= N
            % 如果周期比帧长还长，退化为 Lag-1 检测 (极少见)
            r = xcorr(band_sig, 1, 'coeff'); 
            corr = r(2);
        else
            vec_curr = band_sig(start_idx:end);
            vec_prev = band_sig(1:end-pitch_lag);
            
            numerator = sum(vec_curr .* vec_prev);
            denom = sqrt(sum(vec_curr.^2) * sum(vec_prev.^2)) + 1e-10;
            
            corr = numerator / denom;
        end
        
        % C. 阈值判定 (混合音的关键)
        % 低频带 (0-500, 500-1000) 应该非常严格，容易是浊音
        % 高频带 (2000-4000) 如果是浊音，要求相关性非常强
        
        % 标准建议阈值: 0.6
        % 为了获得更丰富的混合状态，可以对高频带稍微放宽一点，或者保持 0.6
        if corr > 0.6
            bpvc(k) = 1;
        else
            bpvc(k) = 0;
        end
    end
    
    % D. 强制最低频带 (0-500Hz)
    % 如果全局被判为浊音 (调用此函数的前提)，最低频带几乎必须是浊音
    % 这能避免出现 0 1 1 1 1 这种不物理的状态
    bpvc(1) = 1;
    
    % E. 平滑约束 (可选)
    % 强行修正“中间断层”，例如 1 0 1 0 0 -> 1 0 0 0 0
    % MELP 标准有时会强制 BPVC 必须是从低到高的连续 1
    % 但为了保留细节，我们先不做硬性截断
end

