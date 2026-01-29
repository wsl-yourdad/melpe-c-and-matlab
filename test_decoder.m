% =========================================================================
% MELPe 1200bps 解码器 (最终配套版)
% =========================================================================
clc; clear; close all;

% 1. 加载资源
fprintf('加载码本和编码数据...\n');
load('THE_FINAL_CODEBOOK.mat'); % codebooks
load('encoder_output_new.mat'); % all_bit_streams, FS

% 2. 初始化
num_superframes = length(all_bit_streams);
total_frames = num_superframes * 3;
decoded_signal = [];

% 状态记忆
mem_lpc = zeros(10, 1);
mem_disp = zeros(120, 1); % 色散滤波器记忆
seed = 1;
prev_lsf = codebooks.msvq_mean(:); % 初始 LSF (Float)
prev_gain = 10.0; % 初始增益

% 色散滤波器系数 (默认)
disp_filter = [-0.0054, -0.0125, -0.0234, -0.0385, -0.0574, -0.0789, ...
               -0.1009, -0.1206, -0.1348, -0.1408, -0.1373, -0.1235, ...
               -0.0993, -0.0649, -0.0210, 0.0305, 0.0877, 0.1481, ...
               0.2090, 0.2673, 0.3201, 0.3647, 0.3989, 0.4209, ...
               0.4294, 0.4243, 0.4061, 0.3756, 0.3340, 0.2831, ...
               0.2253, 0.1636, 0.1010, 0.0406, -0.0150, -0.0632]'; 
if isfield(codebooks, 'disp_cof'), disp_filter = codebooks.disp_cof; end

% 调试记录容器
dec_debug = struct('pitch', [], 'gain', [], 'lsf_idx', []);

fprintf('开始解码 %d 个超帧...\n', num_superframes);

for s_idx = 1:num_superframes
    b_str = all_bit_streams{s_idx};
    if length(b_str) ~= 81
        % 丢包处理 (简单重复上一帧，或者静音)
        continue; 
    end
    
    % --- A. 解包 (Unpack) ---
    % 严格对应 pack_bits_c_style: 
    % Sync(1)|UV(3)|Par(1)|Pitch(9)|LSF(42)|Gain(10)|BPVC(6)|FS(8)|Jit(1)
    
    p.sync  = b_str(1);
    p.uv    = b_str(2:4); % String '001' etc.
    p.par   = b_str(5);
    p.pitch = bin2dec(b_str(6:14));
    p.lsf_bits = b_str(15:56); % 42 bits string
    p.gain  = bin2dec(b_str(57:66));
    p.bpvc  = bin2dec(b_str(67:72));
    p.fs    = bin2dec(b_str(73:80));
    p.jit   = b_str(81);
    
    % 记录调试信息
    dec_debug.pitch(s_idx) = p.pitch;
    dec_debug.gain(s_idx)  = p.gain + 1; % 记录为 MATLAB 索引以便对比
    
    % --- B. 解码参数 ---
    
    % 1. Pitch
    if p.pitch == 0
        final_pitch = 0;
    else
        % Log 解码: 10^(1.301 + idx * step)
        % 范围 log10(20)~log10(160) -> 1.30103 ~ 2.20412
        step = (2.20412 - 1.30103) / 512;
        final_pitch = 10^(1.30103 + p.pitch * step);
    end
    
    % 2. Gain (查表)
    % 修正: MATLAB 索引 = bit_val + 1
    g_idx = p.gain + 1;
    if g_idx > size(codebooks.gain, 1), g_idx = size(codebooks.gain, 1); end
    
    % 查表得到 2 个增益 (G1, G2) * 3 帧
    % 码本结构: [G1_f1 G2_f1 G1_f2 G2_f2 G1_f3 G2_f3] (Q8 -> /256.0)
    g_vec = codebooks.gain(g_idx, :); 
    % 转换为线性值: 10^(val/20)
    g_lin = 10.^(g_vec / 20.0);
    
    % 3. LSF (核心修复版)
    [lsf_matrix, lsf_dbg_idx] = decode_lsf_sc1200(p.lsf_bits, p.uv, codebooks, prev_lsf);
    dec_debug.lsf_idx(s_idx) = lsf_dbg_idx; % 记录第一个索引用于对比
    
    % 更新预测器记忆 (取第3帧)
    prev_lsf = lsf_matrix(:, 3);
    
    % 4. BPVC, Jitter, FS (简化处理，直接用)
    bpvc_vec = bitget(p.bpvc, 1:5); % 5-band
    if p.jit == '1', jitter = 0.25; else, jitter = 0.0; end
    
    if p.fs > 0 && isfield(codebooks, 'fsvq_cb')
        fs_idx = p.fs + 1;
        if fs_idx <= size(codebooks.fsvq_cb, 2)
            fsmag = codebooks.fsvq_cb(:, fs_idx);
        else
            fsmag = ones(10,1);
        end
    else
        fsmag = ones(10,1);
    end
    
    % --- C. 逐帧合成 ---
    for f = 1:3
        % 准备当前帧参数
        curr_lsf = lsf_matrix(:, f);
        
        % 简单的 Pitch 插值 (如果需要)
        curr_pitch = final_pitch; 
        
        % Gain: 取 G2 (Speech Gain) 
        % 码本顺序: [G1 F1, G2 F1, G1 F2, G2 F2, G1 F3, G2 F3]
        % 索引: (f-1)*2 + 2
        curr_gain = g_lin( (f-1)*2 + 2 );
        
        % 清音帧强制 Pitch=0
        if curr_gain < 1.0, curr_pitch = 0; end % 简单静音门限
        
        % 调用合成 (稳定版)
        [sig_chk, mem_lpc, mem_disp, seed] = melp_syn_frame_enhanced(...
            curr_lsf, curr_pitch, bpvc_vec, curr_gain, prev_gain, ...
            fsmag, mem_lpc, mem_disp, seed, codebooks, disp_filter);
        
        decoded_signal = [decoded_signal; sig_chk];
        prev_gain = curr_gain;
    end
end

% 保存调试数据供对比脚本使用
save('decode_output.mat', 'dec_debug', 'decoded_signal');
fprintf('✅ 解码完成，已保存 decode_output.mat\n');

% 播放
soundsc(decoded_signal, 8000);
audiowrite('out_decoded_final.wav', decoded_signal/max(abs(decoded_signal)), 8000);

% =========================================================================
% 子函数 1: LSF 解码 (Float 域)
% =========================================================================
function [lsf_out, dbg_idx] = decode_lsf_sc1200(bits, uv_str, cb, prev_lsf)
    lsf_out = zeros(10, 3);
    mean_lsf = cb.msvq_mean(:); % Float
    
    if strcmp(uv_str, '001')
        % --- Mode 1 ---
        idx_base = bin2dec(bits(1:9)) + 1;
        idx_int  = bin2dec(bits(10:13)) + 1;
        % MSVQ indices (4 stages)
        idx_s1   = bin2dec(bits(14:21)) + 1;
        idx_s2   = bin2dec(bits(22:27)) + 1;
        idx_s3   = bin2dec(bits(28:33)) + 1;
        idx_s4   = bin2dec(bits(34:39)) + 1;
        
        dbg_idx = idx_base; % 调试用
        
        % 1. Recover Base
        lsf_base = cb.lsp_uv_9(:, idx_base);
        
        % 2. Recover Residual (4-stage MSVQ)
        res = cb.res_cb{1}(:, idx_s1) + ...
              cb.res_cb{2}(:, idx_s2) + ...
              cb.res_cb{3}(:, idx_s3) + ...
              cb.res_cb{4}(:, idx_s4);
              
        % 3. Prediction
        w = cb.inpCoef(:, idx_int); % 20 dim
        w_f2 = w(1:10);
        w_f3 = w(11:20);
        
        % Frame 2
        pred_f2 = w_f2 .* lsf_base + (1 - w_f2) .* prev_lsf;
        % Frame 3 (combined with residual)
        pred_f3 = w_f3 .* lsf_base + (1 - w_f3) .* prev_lsf;
        lsf_f3 = pred_f3 + res(11:20);
        
        % 插值 Frame 1 (简单线性)
        lsf_f1 = 0.5 * prev_lsf + 0.5 * pred_f2; 
        lsf_f2 = pred_f2 + res(1:10);
        
        lsf_out(:,1) = lsf_f1;
        lsf_out(:,2) = lsf_f2;
        lsf_out(:,3) = lsf_f3;
        
    else
        % --- Mode 0 ---
        % F1
        f1_s1 = bin2dec(bits(1:8)) + 1;
        f1_s2 = bin2dec(bits(9:14)) + 1;
        % F2
        f2_s1 = bin2dec(bits(15:22)) + 1;
        f2_s2 = bin2dec(bits(23:28)) + 1;
        % F3
        f3_s1 = bin2dec(bits(29:36)) + 1;
        f3_s2 = bin2dec(bits(37:42)) + 1;
        
        dbg_idx = f1_s1; % 调试用
        
        % Recover (Mean + VQ1 + VQ2)
        lsf_out(:,1) = mean_lsf + cb.lsp_v_cb{1}(:, f1_s1) + cb.lsp_v_cb{2}(:, f1_s2);
        lsf_out(:,2) = mean_lsf + cb.lsp_v_cb{1}(:, f2_s1) + cb.lsp_v_cb{2}(:, f2_s2);
        lsf_out(:,3) = mean_lsf + cb.lsp_v_cb{1}(:, f3_s1) + cb.lsp_v_cb{2}(:, f3_s2);
    end
end

% =========================================================================
% 子函数 2: 合成 (C标准稳定版 - 无蟋蟀声/无溢出)
% =========================================================================
function [sig, mem_lpc, mem_disp, seed] = melp_syn_frame_enhanced(lsf_q15, pitch, bpvc, gain_curr, gain_prev, fsmag, mem_lpc, mem_disp, seed, cb, disp_filter)
    
    FRAME = 180;
    
    % --- 1. LSF 恢复与 C 标准钳位 ---
    lsf_norm = lsf_q15(:);
    
    % 量级归一化 (防止 Q15 误入)
    if mean(abs(lsf_norm)) > 10 
        lsf_norm = lsf_norm / 32768.0; 
    end
    
    lsf_norm = sort(lsf_norm);
    
    % 最小间隔 (50Hz)
    DELTA = 0.0125;
    for i = 1:9
        if lsf_norm(i+1) - lsf_norm(i) < DELTA
            avg = (lsf_norm(i+1) + lsf_norm(i)) / 2;
            lsf_norm(i)   = avg - DELTA/2;
            lsf_norm(i+1) = avg + DELTA/2;
        end
    end
    
    % 边界钳位
    lsf_norm(lsf_norm < DELTA) = DELTA;
    lsf_norm(lsf_norm > (1.0 - DELTA)) = 1.0 - DELTA;
    lsf_norm = sort(lsf_norm); 
    
    % --- 2. 转换为 LPC ---
    lpc_a = lsf2poly(lsf_norm * pi); 
    lpc_a = lpc_a(:);
    
    % --- 3. 带宽扩展 (0.994) - 核心稳健性 ---
    GAMMA = 0.994;
    bw_factor = GAMMA .^ (0:length(lpc_a)-1);
    lpc_a = lpc_a .* bw_factor(:);
    
    % --- 4. 激励生成 ---
    pulse_exc = zeros(1, FRAME);
    if pitch > 0
        T = round(pitch);
        if T < 20, T = 20; end
        pulse_exc(1:T:FRAME) = 1;
        pulse_exc = pulse_exc * sqrt(T);
        mem_disp = mem_disp(:); disp_filter = disp_filter(:);
        [pulse_out, mem_disp] = filter(disp_filter, 1, pulse_exc(:), mem_disp);
        pulse_exc = pulse_out(:).'; 
    else
        mem_disp = zeros(size(mem_disp));
    end
    
    rng(seed);
    noise_exc = randn(1, FRAME);
    seed = rand;
    
    if pitch > 0
        F_pulse = fft(pulse_exc, FRAME); 
        F_noise = fft(noise_exc, FRAME);
        mask = zeros(1, FRAME);
        edges_Hz = [0, 500, 1000, 2000, 3000, 4000];
        edges_bin = round(edges_Hz / (8000/FRAME)); 
        for b = 1:5
            if bpvc(b) == 1
                idx_s = max(2, edges_bin(b)+1); idx_e = edges_bin(b+1);
                if idx_e >= idx_s
                    mask(idx_s : idx_e) = 1;
                    mask(FRAME - idx_e + 2 : FRAME - idx_s + 2) = 1;
                end
            end
        end
        if bpvc(1), mask(1)=1; end
        F_mixed = F_pulse .* mask + F_noise .* (1 - mask);
        excitation = real(ifft(F_mixed));
    else
        excitation = noise_exc;
    end
    
    % --- 5. 增益应用 ---
    g_start_dB = 20*log10(gain_prev + 1e-9);
    g_end_dB   = 20*log10(gain_curr + 1e-9);
    g_traj_dB  = linspace(g_start_dB, g_end_dB, FRAME);
    g_traj     = 10.^(g_traj_dB / 20);
    excitation = excitation .* g_traj;
    
    % --- 6. 合成滤波 (直接滤波) ---
    [sig, mem_lpc] = filter(1, lpc_a, excitation, mem_lpc);
end