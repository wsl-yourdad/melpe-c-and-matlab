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
wav_file = 'OSR_us_000_0010_8k.wav'; % 建议使用你的文件名

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


%  六。主要部分：核心编码循环 (已集成所有C标准模块)
% =========================================================================

% --- 初始化C标准编码结果的容器 ---
fprintf('--- 初始化C标准编码结果容器 ---\n');
superframe_lsf_indices_cell = cell(1, num_superframes);
superframe_pitch_idx_c      = zeros(1, num_superframes);
superframe_bpvc_idx_c       = zeros(1, num_superframes);
superframe_jitter_idx_c     = zeros(1, num_superframes);
superframe_gain_indices     = zeros(1, num_superframes); % 保留
all_bit_streams             = cell(1, num_superframes);

% --- 用于收集子帧数据的临时缓冲区 ---
lpc_a_superframe = zeros(LPC_ORD, 3);

% =========================================================================

% --- 初始化 ---
fprintf('--- 初始化编码器状态与C标准容器 ---\n');
last_lsf_quant = codebooks.msvq_mean(:) / (FS/2); 

num_total_frames = num_superframes * 3;
viz_pitches = zeros(1, num_total_frames);
viz_voicing = zeros(1, num_total_frames);
viz_gains   = zeros(1, num_total_frames);
viz_fsmags  = zeros(10, num_total_frames);
viz_p1      = zeros(1, num_total_frames);
viz_p2      = zeros(1, num_total_frames);
viz_g1      = zeros(1, num_total_frames);

superframe_lsf_indices_cell = cell(1, num_superframes);
superframe_pitch_idx_c      = zeros(1, num_superframes);
superframe_bpvc_idx_c       = zeros(1, num_superframes);
superframe_jitter_idx_c     = zeros(1, num_superframes);
superframe_gain_indices     = zeros(1, num_superframes);
all_bit_streams             = cell(1, num_superframes);

% ====================[ 主循环开始 ]====================
for s_idx = 1:num_superframes
    
    sf_lsfs         = zeros(LPC_ORD, 3);
    sf_pitches      = zeros(1, 3);
    sf_v_dec        = zeros(1, 3);
    sf_gains        = zeros(2, 3);
    sf_fsmags_tmp   = zeros(10, 3);
    lpc_a_superframe = zeros(LPC_ORD, 3);
    
    for f_in_s = 1:3
        frame_idx = (s_idx-1) * 3 + f_in_s;
        start_pt = (frame_idx-1) * FRAME_LEN + 1;
        curr_frame = signal_final(start_pt : start_pt + FRAME_LEN - 1);

        [a, res_frame] = melp_lpc_analysis(curr_frame, LPC_ORD);
        lpc_a_superframe(:, f_in_s) = a;
        sf_lsfs(:, f_in_s) = lpc_to_lsf(a);

        [P1, ~] = melp_pitch_integer(curr_frame, FS);
        [P2, Vp2] = melp_pitch_fraction(curr_frame, P1, FS);
        [P3, ~, Voicing_flag] = melp_pitch_final(curr_frame, FS, P2, Vp2);
        sf_pitches(f_in_s) = P3;
        sf_v_dec(f_in_s) = Voicing_flag;

        [G1, G2] = melp_gain_calculator(curr_frame*32767, res_frame*32767, P3, Voicing_flag);
        sf_gains(:, f_in_s) = [G1; G2];

        if Voicing_flag == 1
            [fsmag_norm, ~] = melp_harmonic_magnitudes(a, P3, FS);
            sf_fsmags_tmp(:, f_in_s) = fsmag_norm;
        else
            sf_fsmags_tmp(:, f_in_s) = ones(10, 1);
        end
        
        global_f_idx = (s_idx-1)*3 + f_in_s;
        viz_p1(global_f_idx) = P1;
        viz_p2(global_f_idx) = P2;
        viz_g1(global_f_idx) = G1;
        viz_pitches(global_f_idx) = P3;
        viz_voicing(global_f_idx) = Voicing_flag;
        viz_gains(global_f_idx) = G2;
        viz_fsmags(:, global_f_idx) = sf_fsmags_tmp(:, f_in_s);
    end
    
    % ========== [超帧级别处理] ==========

    lpc_coeffs_full = [ones(1, 3); lpc_a_superframe];
    [p_idx_c, bpvc_idx_c, jit_idx_c] = quantize_pitch_bpvc_c_style(sf_pitches, sf_v_dec, lpc_coeffs_full, codebooks);
    superframe_pitch_idx_c(s_idx) = p_idx_c;
    superframe_bpvc_idx_c(s_idx) = bpvc_idx_c;
    superframe_jitter_idx_c(s_idx) = jit_idx_c;

    target_gain_vec_6D = reshape(sf_gains, 6, 1);
    joint_gain_idx = gain_vq_search(target_gain_vec_6D, codebooks.gain);
    superframe_gain_indices(s_idx) = joint_gain_idx;
    
    lsf_indices_c = quantize_lsf_c_style(sf_lsfs, sf_v_dec, codebooks, last_lsf_quant);
    superframe_lsf_indices_cell{s_idx} = lsf_indices_c;
    
    fourier_vq_idx = 1;
    last_voiced_frame_idx = find(sf_v_dec == 1, 1, 'last');
    if ~isempty(last_voiced_frame_idx)
         fsmag_target_vec = sf_fsmags_tmp(:, last_voiced_frame_idx);
         if isfield(codebooks, 'fsvq_cb')
            fourier_vq_idx = melp_vq_pitch_search(fsmag_target_vec, codebooks.fsvq_cb);
         end
    end
    
    all_bit_streams{s_idx} = pack_bits_c_style(...
        superframe_lsf_indices_cell{s_idx}, ...
        sf_v_dec, ... % <-- 修正：传入sf_v_dec
        superframe_pitch_idx_c(s_idx), ...
        superframe_bpvc_idx_c(s_idx), ...
        superframe_jitter_idx_c(s_idx), ...
        superframe_gain_indices(s_idx), ...
        fourier_vq_idx - 1 ...
    );
    
end
% ====================[ 主循环结束 ]====================
 %% =========================================================================***********************************************************************************************
%核心循环第二板块： Step 2: Pitch, BPVC, Jitter, and Gain Quantization (C Standard)
% =========================================================================

% --- 2A. Pitch, BPVC, Jitter Quantization (New C-Standard Logic) ---
% 准备LPC系数用于BPVC分析
% 在前面的循环中，我们已经通过 melp_lpc_analysis 得到了系数 'a'
% 我们需要将它们收集起来
%if f_in_s ==
%1^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%    lpc_a_superframe = zeros(LPC_ORD, 3);
%end
%lpc_a_superframe(:, f_in_s) = a; % 'a' 是 melp_lpc_analysis 的输出，不含a0
% -- 在超帧的最后一个子帧循环结束后，执行联合量化
% --^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   if f_in_s == 3
    % 构造完整的LPC系数，包含 a0=1
    lpc_coeffs_full = [ones(1, 3); lpc_a_superframe];

    [p_idx_c, bpvc_idx_c, jit_idx_c] = quantize_pitch_bpvc_c_style(sf_pitches, sf_v_dec, lpc_coeffs_full, codebooks);
    
    superframe_pitch_idx_c(s_idx) = p_idx_c;
    superframe_bpvc_idx_c(s_idx) = bpvc_idx_c;
    superframe_jitter_idx_c(s_idx) = jit_idx_c;


    % --- 2B. Gain Quantization (Unaltered, as requested) ---
    % 增益量化部分保持不变，因为您确认它是正确的 (10-bit VQ)
    target_gain_vec_6D = reshape(sf_gains, 6, 1);
    joint_gain_idx = gain_vq_search(target_gain_vec_6D, codebooks.gain);
    superframe_gain_indices(s_idx) = joint_gain_idx;
   end

%% =========================================================================
% [新] 任务四: 最终比特流打包 (C语言标准)
% =========================================================================
% 同样，在超帧的最后一个子帧循环结束后执行

    if f_in_s == 3
    % 从fourier量化中获取索引
    fourier_vq_idx = 1; % 默认值
    last_voiced_frame_idx = find(sf_v_dec == 1, 1, 'last');
       if ~isempty(last_voiced_frame_idx)
         fsmag_target_vec = sf_fsmags_tmp(:, last_voiced_frame_idx);
         % 假设fsvq_cb码本存在
            if isfield(codebooks, 'fsvq_cb')
            fourier_vq_idx = melp_vq_pitch_search(fsmag_target_vec, codebooks.fsvq_cb);
            end
       end

    all_bit_streams{s_idx} = pack_bits_c_style(...
        superframe_lsf_indices_cell{s_idx}, ...
        sf_v_dec, ... % <-- 修正：传入sf_v_dec
        superframe_pitch_idx_c(s_idx), ...
        superframe_bpvc_idx_c(s_idx), ...
        superframe_jitter_idx_c(s_idx), ...
        superframe_gain_indices(s_idx), ... % 旧的增益索引
        fourier_vq_idx - 1 ... % FS索引, 转为0-based
    );
   end



    %% 4. 调试打印 (抽检)
    %if mod(s_idx, 3) == 0
    %    fprintf('超帧 %2d | 模式: %5s | 联合基音索引: %d | 增益原值: [%.1f, %.1f, %.1f] | 联合增益索引: %d\n', ...
    %            s_idx, pitch_mode, joint_pitch_idx, target_gain_vec_6D(1), target_gain_vec_6D(2), target_gain_vec_6D(3), joint_gain_idx);
    %end

    %fprintf('超帧 %2d | 模式: %5s | 联合基音索引: %d | 增益原值: [%.1f, %.1f, %.1f] | 联合增益索引: %d\\n', ...
    %    s_idx, pitch_mode, joint_pitch_idx, target_gain_vec_6D(1), target_gain_vec_6D(2), target_gain_vec_6D(3), joint_gain_idx);

    %%%%%%%%%%%%%%%加入傅里叶级数幅值矢量量化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % =========================================================================
% 新增模块：傅里叶级数幅值 (Fourier Magnitudes) 矢量量化 (8-bit)
% =========================================================================
fourier_vq_idx = 1; % 默认索引为1 (对应0-based的索引0)

% 1. 寻找超帧中的最后一个浊音帧
last_voiced_frame_idx = find(sf_v_dec == 1, 1, 'last');

   if ~isempty(last_voiced_frame_idx)
    % 2. 如果存在浊音帧, 则使用该帧的傅里叶幅值作为量化目标
    fsmag_target_vec = sf_fsmags_tmp(:, last_voiced_frame_idx);
    
    % 3. 执行VQ搜索, codebooks.fsvq_cb 是 10x256 的码本
    fourier_vq_idx = melp_vq_pitch_search(fsmag_target_vec, codebooks.fsvq_cb);
   else
    % 4. 如果整个超帧都是清音, 索引通常设为0
    fourier_vq_idx = 1; % MATLAB的索引是1-based
   end

%% =========================================================================&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Step 3: LSF 量化 (采用C语言标准的42-bit MSVQ)
% 调用新的多模式量化引擎
lsf_indices_c = quantize_lsf_c_style(sf_lsfs, sf_v_dec, codebooks, last_lsf_quant);

% 将结果存入一个cell，因为它的长度是动态的
superframe_lsf_indices_cell{s_idx} = lsf_indices_c;
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




% ========== [超帧级别处理 - 唯一的正确版本] ==========

% --- Step 2: Pitch, BPVC, Jitter, Gain Quantization ---
lpc_coeffs_full = [ones(1, 3); lpc_a_superframe];
[p_idx_c, bpvc_idx_c, jit_idx_c] = quantize_pitch_bpvc_c_style(sf_pitches, sf_v_dec, lpc_coeffs_full, codebooks);
superframe_pitch_idx_c(s_idx) = p_idx_c;
superframe_bpvc_idx_c(s_idx) = bpvc_idx_c;
superframe_jitter_idx_c(s_idx) = jit_idx_c;

target_gain_vec_6D = reshape(sf_gains, 6, 1);
joint_gain_idx = gain_vq_search(target_gain_vec_6D, codebooks.gain);
superframe_gain_indices(s_idx) = joint_gain_idx;

% --- Step 3: LSF Multi-Mode Quantization ---
lsf_indices_c = quantize_lsf_c_style(sf_lsfs, sf_v_dec, codebooks, last_lsf_quant);
superframe_lsf_indices_cell{s_idx} = lsf_indices_c;

% --- Step 4: Fourier Series Quantization ---
fourier_vq_idx = 1;
last_voiced_frame_idx = find(sf_v_dec == 1, 1, 'last');
if ~isempty(last_voiced_frame_idx)
     fsmag_target_vec = sf_fsmags_tmp(:, last_voiced_frame_idx);
     if isfield(codebooks, 'fsvq_cb')
        fourier_vq_idx = melp_vq_pitch_search(fsmag_target_vec, codebooks.fsvq_cb);
     end
end

% --- Step 5: Final Bit Packing (C-Standard) ---
all_bit_streams{s_idx} = pack_bits_c_style(...
    superframe_lsf_indices_cell{s_idx}, ...
    sf_v_dec, ...
    superframe_pitch_idx_c(s_idx), ...
    superframe_bpvc_idx_c(s_idx), ...
    superframe_jitter_idx_c(s_idx), ...
    superframe_gain_indices(s_idx), ...
    fourier_vq_idx - 1 ...
);

    % 存储到全局比特矩阵 (用于计算总码率)
    full_frame_bits = all_bit_streams{s_idx};


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
    if s_idx <= length(superframe_lsf_indices_cell) && ~isempty(superframe_lsf_indices_cell{s_idx})
        temp_indices = superframe_lsf_indices_cell{s_idx};
        lsf_indices_str = sprintf('%d,%d..', temp_indices(1), temp_indices(2));
    end

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
% =========================================================================
% 9. [最终] 保存编码器输出以供解码器使用
% =========================================================================
fprintf('\n正在为解码器保存输出...\n');
save('encoder_output.mat', 'all_bit_streams', 'signal_final', 'FS');
fprintf('✅ encoder_output.mat 已保存。\n');
