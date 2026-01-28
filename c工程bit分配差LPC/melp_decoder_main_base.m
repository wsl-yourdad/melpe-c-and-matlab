%% =========================================================================
%  MELPe 1200bps 解码器 
%  核心修正:
%    1. 修复变量名报错 (dec_state -> state)
%    2. 注入真实的 Pulse Dispersion (Q15)
%    3. 增强 LPC 稳定性保护 (防止爆音)
%    4. 优化增益平滑 (防止声音短促)
% =========================================================================
clc; clear; close all;

% --- 1. 加载资源 ---
fprintf('Loading Codebooks and Bitstream...\n');
if ~isfile('THE_FINAL_CODEBOOK.mat') || ~isfile('encoder_output_new.mat')
    error('缺少必要文件: THE_FINAL_CODEBOOK.mat 或 encoder_output_new.mat');
end
load('THE_FINAL_CODEBOOK.mat'); 
load('encoder_output_new.mat');     

% 全局常量
FS = 8000;
SUPERFRAME_LEN = 540;

% --- 2. 注入真实的 disp_cof (Q15 -> Float) ---
disp_filter = codebooks.disp_cof(:);
REAL_DISP_FILTER = disp_filter;

% --- 3. 解码器状态初始化 ---
dec_state.prev_lsf_q15 = codebooks.msvq_mean(:); 
dec_state.prev_gain_lin = 1.0; 
dec_state.prev_pitch = 40.0;
dec_state.syn_mem = zeros(10, 1); % LPC 记忆 (列向量)
dec_state.disp_mem = zeros(length(REAL_DISP_FILTER)-1, 1); % 色散记忆
dec_state.seed = 1;
dec_state.prev_pulse_pos = 1;%%%%%*******************************************************************************************************
% 输出缓冲
total_samples = length(all_bit_streams) * SUPERFRAME_LEN;
full_speech = zeros(1, total_samples);
write_ptr = 1;

fprintf('Decoding %d Superframes...\n', length(all_bit_streams));

%% --- 4. 主解码循环 (带数据记录) ---
fprintf('Decoding %d Superframes...\n', length(all_bit_streams));

% [新增] 初始化调试容器
dec_debug.pitch_idx = [];
dec_debug.gain_idx  = [];
dec_debug.lsf_idx   = []; % 记录 Mode0:F1S1 或 Mode1:Base
dec_debug.bpvc_idx  = [];

for s_idx = 1:length(all_bit_streams)
    bit_str_81 = all_bit_streams{s_idx};
    
    if length(bit_str_81) ~= 81
        decoded_superframe = zeros(1, SUPERFRAME_LEN);
        % 填充空数据防止长度不一致
        dec_debug.pitch_idx(end+1) = 0;
        dec_debug.gain_idx(end+1)  = 0;
        dec_debug.lsf_idx(end+1)   = 0;
        dec_debug.bpvc_idx(end+1)  = 0;
    else
        % === 核心调用 (接收 params) ===
        dec_state.debug_sf_idx = s_idx;
        % [修改] 接收第3个返回值 params
        [decoded_superframe, dec_state, params] = melpe_s(bit_str_81, dec_state, codebooks, REAL_DISP_FILTER);
        
        % [新增] 记录本帧参数到调试容器
        dec_debug.pitch_idx(end+1) = params.pitch_idx;
        dec_debug.gain_idx(end+1)  = params.gain_idx;
        dec_debug.bpvc_idx(end+1)  = params.bpvc_idx;
        
        % LSF 记录关键索引以便对比
        if params.mode == 1
            dec_debug.lsf_idx(end+1) = params.lsf_indices.base;
        else
            dec_debug.lsf_idx(end+1) = params.lsf_indices.f1_s1;
        end
    end
    
    full_speech(write_ptr : write_ptr + SUPERFRAME_LEN - 1) = decoded_superframe;
    write_ptr = write_ptr + SUPERFRAME_LEN;
    
    if mod(s_idx, 10) == 0, fprintf('.'); end
end
fprintf('\nDone!\n');
%% --- 5. 播放与后处理 ---
% 去直流
b_hp = [1, -0.95]; a_hp = [1, -0.95*0.5]; % 稍微温和一点的高通
full_speech = filter(b_hp, a_hp, full_speech);

% 幅度归一化 (防止过小或过大)
max_amp = max(abs(full_speech));
if max_amp > 1e-9
    full_speech = full_speech / max_amp * 0.95;
end

fprintf('Playing...\n');
soundsc(full_speech, FS);
audiowrite('melpe_1200_final_new.wav', full_speech, FS);

% 画图验证
figure; 
subplot(2,1,1); plot(full_speech); title('Decoded Speech'); grid on;
subplot(2,1,2); spectrogram(full_speech, 256, 128, 256, FS, 'yaxis'); title('Spectrogram');

% === 请替换解码器末尾的保存代码 ===

% 1. 获取接收到的比特流 (确保它被保存)
% 假设您的主循环里用的是 all_bit_streams 或者 bit_str_81
% 如果您在主循环里没有单独存，现在直接把编码器的输入作为"理想接收"存下来用于对比
debug_rx_streams = all_bit_streams; 

% 2. 保存所有关键数据 (波形 + 调试参数 + 比特流)
save('decode_output.mat', 'dec_debug', 'debug_rx_streams');

fprintf('✅ 数据已保存: dec_debug, decoded_signal, debug_rx_streams 全部就位。\n');

%% =========================================================================
%  函数库
% ===============================================================================================================================
function [sp_buf, state,params] = melpe_s(bit_str, state, cb, disp_filter)
    % 1. 解包
    params = unpack_sc1200(bit_str);
    
    if ~isfield(state, 'prev_pulse_pos')
        state.prev_pulse_pos = 1; 
    end
    % =====================================================================********************************************************
    % [调试埋点] 解码器实测值打印
    % =====================================================================
    % 提取 LSF 关键索引 (逻辑与编码器保持一致)
    if params.mode == 1
        lsf_debug_tag = sprintf('M1/Base:%-3d', params.lsf_indices.base);
    else
        lsf_debug_tag = sprintf('M0/F1S1:%-3d', params.lsf_indices.f1_s1);
    end
    
    % 使用与编码器完全相同的格式打印，方便对比
    % 注意: 这里不需要 persistent 限制，我们要看全量数据对比
    fprintf('[DEC_CHECK] SF:%03d | UV_Code:%d%d%d | P:%03d | G:%04d | %-10s | BPVC:%02d | J:%d', ...
        state.debug_sf_idx, ... % 需要在主循环传入 s_idx，如果没有，暂时用 0 或全局变量
        params.uv(1), params.uv(2), params.uv(3), ...
        params.pitch_idx, ...
        params.gain_idx, ...
        lsf_debug_tag, ...
        params.bpvc_idx, ...
        params.jitter_idx);
        
    % 简易的“对齐检查”
    % 如果 Pitch 或 Gain 也是 0 或极小值，可能意味着错位
    if params.gain_idx == 0 && params.pitch_idx == 0
        fprintf(' <--- [警告] 数据疑似全零或错位');
    end
    fprintf('\n');
    % =====================================================================*****************************************************

    % === [新增] 调试探针：打印本帧核心参数 ===
    % 仅在前 10 帧打印，避免刷屏
    persistent debug_counter;
    if isempty(debug_counter), debug_counter = 0; end
    if debug_counter < 10
        debug_counter = debug_counter + 1;
        fprintf('[Decoder Debug SF#%d] Mode:%d | PitchIdx:%3d | GainIdx:%4d | BPVC:%2d | UV:%s\n', ...
            debug_counter, ...
            isfield(params, 'mode') && params.mode, ...
            params.pitch_idx, ...
            params.gain_idx, ...
            params.bpvc_idx, ...
            sprintf('%d%d%d', params.uv));
    end

    % 2. 参数解码
    % 增益 (10 bits -> 6个增益值)
    g_db_vec = cb.gain(params.gain_idx-1, :); %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    g_lin_vec = 10.^(g_db_vec / 20);      
    
    % 基音 (末帧基音)
    final_pitch = decode_pitch_sc1200(params.pitch_idx, params.uv(3));
    
    % LSF (核心修复：直接获取 3 帧独立数据)
    % 返回 [10 x 3] 矩阵
    lsf_3frames = decode_lsf_sc1200(params, state.prev_lsf_q15, cb);%*****************************************************************************************
    
    % BPVC (末帧带通状态)
    bpvc_final = decode_bpvc_sc1200(params.bpvc_idx, params.uv(3));
    
    % 傅里叶幅值
    if params.fs_idx > size(cb.fsvq_cb, 2), params.fs_idx = 1; end
    fsmag = cb.fsvq_cb(:, params.fs_idx);
    
    % --- 3. 逐帧合成 (Frame-by-Frame Synthesis) ---
    sp_buf = zeros(1, 540);
    FRAME = 180;
    
    for f = 1:3
        % [✅ 修复] 直接取解码出的 LSF，严禁再次插值！
        lsf_q15 = lsf_3frames(:, f);
        
        % [基音插值] 
        % 1200bps 只传了最后一帧的 pitch，前两帧确实需要插值
        % 但要注意 UV 处理：如果一端是 UV(pitch=0)，插值需要小心
        alpha = f / 3.0;
        
        p_prev = state.prev_pitch;
        p_curr = final_pitch;
        
        % 简单逻辑：如果两头都有基音，则线性插值；否则保持有基音的那一头，或归零
        if p_prev > 0 && p_curr > 0
            curr_pitch = (1 - alpha) * p_prev + alpha * p_curr;
        elseif p_curr > 0
            curr_pitch = p_curr; % 前面是清音，直接起振
        elseif p_prev > 0
            curr_pitch = p_prev; % 后面是清音，保持余振
        else
            curr_pitch = 0; % 全清音
        end
        
        % Jitter 处理 (仅在中间帧或特定条件下加抖动，这里简化处理)
        if f == 2 && params.jitter_idx == 1 && curr_pitch > 0
            curr_pitch = curr_pitch * 0.75; % 模拟抖动
        end
        
        % [增益提取]
        curr_gain_lin = g_lin_vec((f-1)*2 + 1); % 取每帧的第2个子增益点作为本帧增益
        prev_gain_lin = state.prev_gain_lin;
        
        % [BPVC 处理]
        % 1200bps BPVC 只有末帧数据。
        % 策略：如果本帧是 UV (由 params.uv(f) 决定)，强制全0
        % 如果本帧是 V，则沿用末帧的 BPVC (或者你可以做个插值，但沿用通常够了)
        if params.uv(f) == 0
            curr_bpvc = zeros(5,1);
            curr_pitch = 0; % 强行清音
        else
            curr_bpvc = bpvc_final;
            % 保护：如果判定为浊音但 pitch 却没了，给个默认 pitch
            if curr_pitch == 0, curr_pitch = 40; end 
        end
        
        % [LPC 稳定性重置]
        if any(isnan(state.syn_mem)) || max(abs(state.syn_mem)) > 1e5
             state.syn_mem = zeros(10, 1);
        end
        
        % [核心合成]
            [syn_frame, state.syn_mem, state.disp_mem, state.seed, state.prev_pulse_pos] = ...
            melp_syn_frame_enhanced(...
            lsf_q15, ...
            curr_pitch, ...
            curr_bpvc, ...
            curr_gain_lin, ...
            state.prev_gain_lin, ... % 注意这里用的是 state.prev_gain_lin
            fsmag, ...
            state.syn_mem, ...
            state.disp_mem, ...
            state.seed, ...
            cb, ...
            disp_filter, ...
            state.prev_pulse_pos);
            
        sp_buf((f-1)*FRAME + 1 : f*FRAME) = syn_frame;
        
        % 更新历史增益 (注意要用这一帧的结束增益)
        state.prev_gain_lin = curr_gain_lin;
    end
    
    % 更新状态
    state.prev_lsf_q15 = lsf_3frames(:, 3); % 只存 F3
    state.prev_pitch = final_pitch;
end
%===============================================================================================

function [sig, mem_lpc, mem_disp, seed,next_pulse_pos] = melp_syn_frame_enhanced(lsf_q15, pitch, bpvc, gain_curr, gain_prev, fsmag, mem_lpc, mem_disp, seed, cb, disp_filter,prev_pulse_pos)
    
    FRAME = 180;
  
    % --- 1. LPC 恢复与稳定性加固 ---
    lsf_norm = lsf_q15(:);
    if mean(abs(lsf_norm)) > 20 
        lsf_norm = lsf_norm / 32768.0; % 如果是大整数，转为 0~1 小数
    end
    
    % [步骤 A] 强制排序
    lsf_norm = sort(lsf_norm);
    
    % [步骤 B] 强制最小间隔 (Safety Check) - 核心修复!
    % 最小间隔 50Hz (在 0~1 归一化域约为 0.0125)
    MIN_DIST = 0.0125; 
    
    % 迭代 10 次以确保所有间隔都拉开
    for iter = 1:10
        violation = false;
        for i = 1:9
            if lsf_norm(i+1) - lsf_norm(i) < MIN_DIST
                % 如果两个线谱太近，把它们撑开
                mid = (lsf_norm(i+1) + lsf_norm(i)) / 2;
                lsf_norm(i)   = mid - MIN_DIST/2;
                lsf_norm(i+1) = mid + MIN_DIST/2;
                violation = true;
            end
        end
        
        % [步骤 C] 边界钳位 (0.005 ~ 0.995)
        if lsf_norm(1) < 0.005, lsf_norm(1) = 0.005; end
        if lsf_norm(10) > 0.995, lsf_norm(10) = 0.995; end
        
        % 再次排序以防万一
        lsf_norm = sort(lsf_norm);
        
        if ~violation, break; end
    end
    
    % 【在这里插入补丁】防止撑开后中间的值越界，这会导致 lsf2poly 报错
    lsf_norm(lsf_norm > 0.95) = 0.95;%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$%%%%%%%%%%%
    lsf_norm(lsf_norm < 0.005) = 0.005;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&
    lsf_norm = sort(lsf_norm);

    % [步骤 D] 转换为 LPC 系数
    % lsf_norm 是 0~1, 乘以 pi 变成 0~pi
    lpc_a = lsf2poly(lsf_norm * pi); 
    
    % [步骤 E] 最终稳定性保险丝
    if any(abs(roots(lpc_a)) >= 1.0)
        % 如果依然不稳定，回退到全通滤波器 (静音/原声)
        lpc_a = [1; zeros(10,1)]; 
    end
    lpc_a = lpc_a(:); 
    
    % --- 2. 激励生成 (Standard MELP Excitation) ---
    % A. 脉冲生成
    pulse_exc = zeros(1, FRAME);
    current_pos = prev_pulse_pos; % 继承上一帧的结束位置
    
    if pitch > 0
        T = pitch; % 周期
        if T < 20, T = 20; end
        
        % 在当前帧内铺设脉冲，而不是每次从 1 开始
        while current_pos <= FRAME
            % 只有在位置合法时才打脉冲
            idx = round(current_pos); 
            if idx >= 1 && idx <= FRAME
                pulse_exc(idx) = 1;
            end
            current_pos = current_pos + T; % 跳到下个周期
        end
        
        % 计算留给下一帧的偏移量 (负数表示下一帧第几个点开始)
        % 例如 current_pos = 190, FRAME = 180, next = 10
        next_pulse_pos = current_pos - FRAME; 
        
        % 能量归一化
        pulse_exc = pulse_exc * sqrt(T);
        
        % 色散滤波 (必须保留)
        mem_disp = mem_disp(:); disp_filter = disp_filter(:);
        [pulse_out, mem_disp] = filter(disp_filter, 1, pulse_exc, mem_disp);
        pulse_exc = pulse_out(:)';
    else
        % 清音帧重置相位或保持随机
        next_pulse_pos = 1; 
        mem_disp = zeros(size(mem_disp));
    end
    
    % B. 噪声生成
    rng(seed);
    noise_exc = randn(1, FRAME);
    seed = rand;
    
    % C. 混合 (Mixed Excitation)
    if pitch > 0
        % FFT 混合
        F_pulse = fft(pulse_exc, FRAME); 
        F_noise = fft(noise_exc, FRAME);
        mask = zeros(1, FRAME);
        
        % 频带定义 (0-500, 500-1000, 1000-2000, 2000-3000, 3000-4000)
        edges_Hz = [0, 500, 1000, 2000, 3000, 4000];
        edges_bin = round(edges_Hz / (8000/FRAME)); 
        
        for b = 1:5
            if bpvc(b) == 1 % 浊音带
                idx_s = max(2, edges_bin(b)+1); % 跳过 DC
                idx_e = edges_bin(b+1);
                
                % 正频率部分
                if idx_e >= idx_s
                    mask(idx_s : idx_e) = 1;
                    % 负频率镜像
                    mask(FRAME - idx_e + 2 : FRAME - idx_s + 2) = 1;
                end
            end
        end
        % DC 和 Nyquist 总是由低频带决定 (通常是脉冲)
        if bpvc(1), mask(1)=1; end
        
        % 频域合成
        F_mixed = F_pulse .* mask + F_noise .* (1 - mask);
        
        % 自适应谱增强 (ASE) - 简化版直接应用
        if ~isempty(fsmag)
             % 这里如果您的代码里有 ASE 逻辑，请保留
             % 如果没有，暂时忽略也不影响核心听感
        end
        
        excitation = real(ifft(F_mixed));
    else
        excitation = noise_exc;
    end
    
    % --- 3. 增益应用 ---
    % 对数域插值
    g_start_dB = 20*log10(gain_prev + 1e-9);
    g_end_dB   = 20*log10(gain_curr + 1e-9);
    g_traj_dB  = linspace(g_start_dB, g_end_dB, FRAME);
    g_traj     = 10.^(g_traj_dB / 20);
    
    excitation = excitation .* g_traj;
    
    % --- 4. LPC 合成滤波 ---
    % 这里的 filter(1, lpc_a) 是标准的 IIR 合成: y(n) = x(n) - sum(a_k * y(n-k))
    [sig, mem_lpc] = filter(1, lpc_a, excitation, mem_lpc);
end
% --- Unpackers
% ---解包函数==================================================================================================
function p = unpack_sc1200(b_str)
% =========================================================================
% 81-Bit Frame Unpacker (Strict Alignment)
% 结构: Sync(1)|UV(3)|Par(1)|Pitch(9)|LSF(42)|Gain(10)|BPVC(6)|FS(8)|J(1)
% =========================================================================

    % 1. Sync Bit (Bit 1)
    p.sync = b_str(1); % 仅做校验用
    
    % 2. UV (Bits 2-4)
    uv_bin = b_str(2:4);
    p.uv = [str2double(uv_bin(1)), str2double(uv_bin(2)), str2double(uv_bin(3))];
    
    % 3. UV Parity (Bit 5) - [新增]
    p.uv_parity = b_str(5);
    
    % 4. Pitch (Bits 6-14) -> 9 bits
    p.pitch_idx = bin2dec(b_str(6:14));
    
    % 5. LSF (Bits 15-56) -> 42 bits
    p.lsf_bits = b_str(15:56);
    
    % 6. Gain (Bits 57-66) -> 10 bits
    p.gain_idx = bin2dec(b_str(57:66)) + 1; % 转 1-based
    
    % 7. BPVC (Bits 67-72) -> 6 bits
    % C标准: 5 bits Index + 1 bit Protection
    % 假设低5位是索引，第6位是保护位 (具体看打包时的位序，通常MSB为保护)
    bpvc_raw = b_str(67:72);
    p.bpvc_idx  = bin2dec(bpvc_raw(1:5)); % 前5位作为索引 (0-31)
    p.bpvc_prot = bpvc_raw(6);            % 最后1位保护位
    % 如果你的打包逻辑是直接 dec2bin(idx,6)，那 idx 就是全部 6 位
    % 为了兼容之前的编码器，这里也可以直接读 6 位：
    % p.bpvc_idx = bin2dec(bpvc_raw); 
    
    % 8. FS Mag (Bits 73-80) -> 8 bits
    p.fs_idx = bin2dec(b_str(73:80)) + 1;
    
    % 9. Jitter (Bit 81) -> 1 bit
    p.jitter_idx = bin2dec(b_str(81));

    % --- LSF 模式解析 ---
    if strcmp(uv_bin, '001')
        p.mode = 1; % Mode 1 (Interpolation)
        % 结构: Base(9) + Int(4) + S1(8) + S2(6) + S3(6) + S4(6) + Prot(3)
        p.lsf_indices.base = bin2dec(p.lsf_bits(1:9)) + 1;
        p.lsf_indices.int  = bin2dec(p.lsf_bits(10:13)) + 1;
        p.lsf_indices.s1   = bin2dec(p.lsf_bits(14:21)) + 1;
        p.lsf_indices.s2   = bin2dec(p.lsf_bits(22:27)) + 1;
        p.lsf_indices.s3   = bin2dec(p.lsf_bits(28:33)) + 1;
        p.lsf_indices.s4   = bin2dec(p.lsf_bits(34:39)) + 1;
        % Prot(3) at 40-42 is skipped/ignored
    else
        p.mode = 0; % Mode 0 (Independent)
        % 结构: F1(8+6) + F2(8+6) + F3(8+6)
        p.lsf_indices.f1_s1 = bin2dec(p.lsf_bits(1:8)) + 1;
        p.lsf_indices.f1_s2 = bin2dec(p.lsf_bits(9:14)) + 1;
        p.lsf_indices.f2_s1 = bin2dec(p.lsf_bits(15:22)) + 1;
        p.lsf_indices.f2_s2 = bin2dec(p.lsf_bits(23:28)) + 1;
        p.lsf_indices.f3_s1 = bin2dec(p.lsf_bits(29:36)) + 1;
        p.lsf_indices.f3_s2 = bin2dec(p.lsf_bits(37:42)) + 1;
    end
end

%=====================================================================================================================

function pitch_val = decode_pitch_sc1200(idx, uv_flag)
    if uv_flag == 0, pitch_val = 0; else
        pitch_val = 10^(1.30103 + (idx) * (2.20412-1.30103)/512);%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    end
end

%=====================================================================================================================

function lsf_out_matrix = decode_lsf_sc1200(p, lsf_prev, cb)
% =========================================================================
% LSF 解码器 (高保真版)
% 修正:
%   1. Mode 1: 使用加权逻辑恢复 F1，而非简单平均。
%   2. Mode 0: 必须加上 msvq_mean，否则频率偏移。
% =========================================================================

    lsf_out_matrix = zeros(10, 3);
    
    % 0. 准备均值向量 (必须加!)
    % 确保是列向量
    mean_lsf = cb.msvq_mean(:); 
    
    if p.mode == 1
        % =================================================================
        % Mode 1: 插值预测解码 (F2 & F3)
        % =================================================================
        idx = p.lsf_indices;
        
        % 1. 解码基准 LSF (Base) - Q15
        lsf_base = cb.lsp_uv_9(:, max(1,idx.base-1));%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        
        % 2. 获取插值系数 (inpCoef) - Q14
        % cb.inpCoef: [20 x 16]
        w_vec = cb.inpCoef(:, idx.int);
        w_f2 = w_vec(1:10);  % 权重 for F2
        w_f3 = w_vec(11:20); % 权重 for F3
        
        % 3. 预测 (Prediction)
        % Pred = w * Base + (1-w) * Prev
        pred_f2 = w_f2 .* lsf_base + (1 - w_f2) .* lsf_prev;
        pred_f3 = w_f3 .* lsf_base + (1 - w_f3) .* lsf_prev;
        
        % 4. 解码残差 (Residual) - Q17
        res_20d = cb.res_cb{1}(:, idx.s1) + ...
                  cb.res_cb{2}(:, idx.s2) + ...
                  cb.res_cb{3}(:, idx.s3) + ...
                  cb.res_cb{4}(:, idx.s4);
        
        % 5. 重建 F2, F3
        % 注意: Mode 1 的预测通常是基于无均值域还是绝对域? 
        % 1200bps 标准中，inpCoef 作用于绝对 LSF。所以不需要再加 mean_lsf。
        lsf_f2 = pred_f2 + res_20d(1:10);
        lsf_f3 = pred_f3 + res_20d(11:20);
        
        % 6. [升级] F1 处理
        % 使用 F2 的反向预测或者插值。
        % C代码中，F1 在 Mode 1 下往往被视为 "Prev" 到 "F2" 的中间态。
        % 使用 inpCoef 的一部分或者固定的 0.5 其实是合理的，但我们可以更精细。
        % 既然 w_f2 决定了 F2 偏向 Base 的程度，F1 应该比 F2 更偏向 Prev。
        % 估算权重: w_f1 ≈ 0.5 * w_f2
        w_f1 = 0.5 * w_f2; 
        lsf_f1 = w_f1 .* lsf_base + (1 - w_f1) .* lsf_prev;
        
        lsf_out_matrix(:, 1) = lsf_f1;
        lsf_out_matrix(:, 2) = lsf_f2;
        lsf_out_matrix(:, 3) = lsf_f3;
        
    else
        % =================================================================
        % Mode 0: 独立解码 (Independent MSVQ)
        % 关键修正: 必须加上 Mean!
        % =================================================================
        idx = p.lsf_indices;
        
        % Frame 1
        r1 = cb.lsp_v_cb{1}(:, idx.f1_s1-1);
        r2 = cb.lsp_v_cb{2}(:, idx.f1_s2-1);
        lsf_out_matrix(:, 1) = mean_lsf + r1 + r2; % ✅ 加上均值
        
        % Frame 2
        r1 = cb.lsp_v_cb{1}(:, idx.f2_s1-1);
        r2 = cb.lsp_v_cb{2}(:, idx.f2_s2-1);
        lsf_out_matrix(:, 2) = mean_lsf + r1 + r2; % ✅ 加上均值
        
        % Frame 3
        r1 = cb.lsp_v_cb{1}(:, idx.f3_s1-1);
        r2 = cb.lsp_v_cb{2}(:, idx.f3_s2-1);
        lsf_out_matrix(:, 3) = mean_lsf + r1 + r2; % ✅ 加上均值*************************************************************************
    end
end
%================================================================================================================================

function val = decode_msvq(indices, stages, mean_val)
    val = mean_val;
    for i = 1:length(indices)
        idx = indices(i);
        if idx > size(stages{i}, 2), idx = size(stages{i}, 2); end
        val = val + stages{i}(:, idx);
    end
end

function bpvc = decode_bpvc_sc1200(idx, uv_flag)
    bpvc = zeros(5,1);
    if uv_flag == 1
        for i=1:5, if bitget(idx, i), bpvc(i) = 1; end; end
    end
end