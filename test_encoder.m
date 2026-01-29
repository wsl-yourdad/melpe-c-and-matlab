%% =========================================================================
%  MELPe 1200bps 编码器 (V15.0 - 脚本对齐查表版)
%  逻辑:
%    1. 直接加载用户生成的 THE_FINAL_CODEBOOK.mat (浮点格式)
%    2. 计算浮点特征 -> 在浮点码本中搜索索引 -> 输出索引
%    3. 绝不自己瞎搞量纲，完全依赖码本的量纲
% =========================================================================
clc; clear; close all;

% --- 1. 加载资源 ---
fprintf('1. Loading Your Codebooks...\n');
if ~isfile('THE_FINAL_CODEBOOK.mat')
    error('请先运行您的 build_final_codebooks.m 脚本生成 .mat 文件！');
end
load('THE_FINAL_CODEBOOK.mat'); % 包含 codebooks 结构体

% 提取码本引用 (方便后续调用)
CB = codebooks;
fprintf('   Codebooks loaded. Mode B stages: %d\n', length(CB.lsf_cb_mode_B));

% --- 2. 音频预处理 ---
[filename, pathname] = uigetfile({'*.wav';'*.mp3'}, 'Select Input Audio');
if isequal(filename,0), error('No file selected'); end
[sig, fs_orig] = audioread(fullfile(pathname, filename));

FS = 8000;
if fs_orig ~= FS, sig = resample(sig, FS, fs_orig); end
if size(sig, 2) > 1, sig = mean(sig, 2); end
% 归一化
sig = sig / max(abs(sig)) * 0.95;

% [C-Style HPF] 这一步必须做，否则 LPC 分析会被低频干扰
b_hp = [1, -1]; a_hp = [1, -0.95];
sig = filter(b_hp, a_hp, sig);
signal_final = sig;

FRAME_LEN = 180;
num_frames = floor(length(sig) / FRAME_LEN);
num_superframes = floor(num_frames / 3);

fprintf('2. Encoding %d Superframes...\n', num_superframes);

% --- 3. 编码主循环 ---
all_bit_streams = {};
LPC_ORD = 10;

% 状态变量 (用于插值)
% 注意：你的脚本里 msvq_mean 是 raw data，但 lsf_cb 是归一化后的吗？
% 看你的脚本：codebooks.msvq_mean = raw_msvq_mean(:); (没有除法)
% 但是 lsf_cb_mode_B 也没有除法。
% 这意味着你的码本里存的都是 Q15 整数！
% 所以我的编码器计算出的 LSF 也必须转成 Q15 整数才能去查表！

prev_lsf_q15 = double(CB.msvq_mean); 

for s_idx = 1:num_superframes
    
    % === Step A: 特征提取 (Hz -> Q15) ===
    sf_lsfs_q15 = zeros(10, 3);
    sf_pitches = zeros(1, 3);
    sf_gains_db = zeros(2, 3);
    sf_uv = zeros(1, 3);
    
    for f = 1:3
        idx_start = (s_idx-1)*3*FRAME_LEN + (f-1)*FRAME_LEN + 1;
        % 加长窗分析 (200 samples matching C)
        ana_len = 200;
        if idx_start + ana_len - 1 <= length(sig)
            sp_ana = sig(idx_start : idx_start + ana_len - 1);
        else
            sp_ana = [sig(idx_start:end); zeros(ana_len-(length(sig)-idx_start+1),1)];
        end
        sp_frm = sp_ana(1:180); 
        
        % 1. Pitch & Voicing
        [P1, ~] = melp_pitch_integer(sp_ana, FS);
        [P2, Vp2] = melp_pitch_fraction(sp_ana, P1, FS);
        [P3, ~, uv_flag] = melp_pitch_final(sp_ana, FS, P2, Vp2);
        sf_pitches(f) = P3;
        sf_uv(f) = uv_flag;
        
        % 2. Gain (dB)
        e_rms = sqrt(mean(sp_frm.^2)) + 1e-9;
        g_db = 20*log10(e_rms);
        sf_gains_db(:, f) = [g_db; g_db];
        
        % 3. LPC -> LSF (Q15 Integer)
        win = hamming(ana_len);
        sp_win = sp_ana .* win;
        r = xcorr(sp_win, LPC_ORD, 'biased');
        r = r(LPC_ORD+1:end);
        a = levinson(r, LPC_ORD);
        
        % [BWEX] 0.994
        bw = 0.994 .^ (0:LPC_ORD);
        a = a .* bw;
        
        lsf_rad = lpc_to_lsf(a, []);
        
        % [Mapping] 0~pi -> 0~32768
        % 你的码本是整数，所以这里必须转整数
        lsf_int = (lsf_rad / pi) * 32768.0;
        
        % [Clamp]
        lsf_int = sort(lsf_int);
        lsf_int(lsf_int < 50) = 50; 
        for k=1:9, if lsf_int(k+1)-lsf_int(k)<50, lsf_int(k+1)=lsf_int(k)+50; end; end
        
        sf_lsfs_q15(:, f) = lsf_int;
    end
    
    % === Step B: 量化 (查表) ===
    
    % 确定模式: 只要最后一帧是 UV，就切 Mode 1 (SC1200 标准)
    anchor_uv = (sf_uv(3) == 0);
    
    bits_lsf = '';
    
    if anchor_uv 
        % === Mode 1 (UV) ===
        % 你的脚本: codebooks.lsf_cb_9bit (512x10)
        % F1, F2: 绝对量化 (9 bits)
        
        % F1
        [~, idx1] = min(sum((double(CB.lsf_cb_9bit)' - sf_lsfs_q15(:,1)).^2, 1));
        % F2
        [~, idx2] = min(sum((double(CB.lsf_cb_9bit)' - sf_lsfs_q15(:,2)).^2, 1));
        
        % F3: MSVQ (24 bits) - 使用 Mode B 的码本 (SC1200通常复用)
        % Target: Residual (LSF - Mean)
        target = sf_lsfs_q15(:,3) - double(CB.msvq_mean);
        [idx_msvq, q_res] = run_msvq(target, CB.lsf_cb_mode_B);
        
        % Update State
        prev_lsf_q15 = double(CB.msvq_mean) + q_res;
        
        % Pack: 9+9+8+6+5+5
        bits_lsf = [dec2bin(idx1-1, 9), dec2bin(idx2-1, 9), ...
                    dec2bin(idx_msvq(1)-1, 8), dec2bin(idx_msvq(2)-1, 6), ...
                    dec2bin(idx_msvq(3)-1, 5), dec2bin(idx_msvq(4)-1, 5)];
        
    else
        % === Mode 2 (V) ===
        % 1. Quantize Anchor F3 (Residual MSVQ)
        target_f3 = sf_lsfs_q15(:,3) - double(CB.msvq_mean);
        % 你的脚本: codebooks.lsf_cb_mode_B (Cell array of 4 matrices)
        [idx_msvq_f3, q_res_f3] = run_msvq(target_f3, CB.lsf_cb_mode_B);
        
        q_f3 = double(CB.msvq_mean) + q_res_f3;
        
        % 2. Interpolate
        pred_f1 = 0.6 * prev_lsf_q15 + 0.4 * q_f3;
        pred_f2 = 0.4 * prev_lsf_q15 + 0.6 * q_f3;
        
        % 3. Quantize Residuals for F1/F2
        res_f1 = sf_lsfs_q15(:,1) - pred_f1;
        res_f2 = sf_lsfs_q15(:,2) - pred_f2;
        
        % F1 Res: 8 bits -> Use Mode B Stage 1 (256)
        cb_s1 = double(CB.lsf_cb_mode_B{1})'; % [10 x 256] -> check transpose in script
        % Script: matrix_1200 is [10 x N], so cb_s1 is [10 x 256]. 
        % Transpose to [10 x 256] for calculation if needed, standard is columns=vectors
        % Let's Assume CB columns are vectors.
        [~, idx_res1] = min(sum((double(CB.lsf_cb_mode_B{1}) - res_f1).^2, 1));
        
        % F2 Res: 6 bits -> Use Mode B Stage 2 (64)
        [~, idx_res2] = min(sum((double(CB.lsf_cb_mode_B{2}) - res_f2).^2, 1));
        
        % Update State
        prev_lsf_q15 = q_f3;
        
        % Pack: 24(F3) + 4(Dummy) + 8(Res1) + 6(Res2)
        bits_lsf = [dec2bin(idx_msvq_f3(1)-1, 8), dec2bin(idx_msvq_f3(2)-1, 6), ...
                    dec2bin(idx_msvq_f3(3)-1, 5), dec2bin(idx_msvq_f3(4)-1, 5), ...
                    '0000', ...
                    dec2bin(idx_res1-1, 8), dec2bin(idx_res2-1, 6)];
    end
    
    % --- Pitch Quantization (9 bits) ---
    % Use CB.pitch_vq_cb_vvv? No, that's VQ (3 frames). SC1200 uses scalar pitch.
    % We use standard Log mapping
    p_val = sf_pitches(3);
    p_log = log10(p_val);
    min_p = 1.3; max_p = 2.3;
    p_idx = round((p_log - min_p)/(max_p - min_p) * 511);
    if p_idx<0, p_idx=0; end; if p_idx>511, p_idx=511; end
    
    % --- Gain Quantization (10 bits) ---
    % Search CB.gain (6xN)
    target_g = sf_gains_db(:)'; % 1x6
    % CB.gain is likely transposed or normalized in your script?
    % Script: codebooks.gain = reshaped_gain_cb' / 256.0;
    % Assuming CB.gain is [N x 6] or [6 x N]. Let's check size
    if size(CB.gain, 1) == 6
        [~, g_idx] = min(sum((CB.gain - target_g').^2, 1));
    else
        [~, g_idx] = min(sum((CB.gain' - target_g').^2, 1));
    end
    
    % --- Assembly ---
    b_uv = [num2str(sf_uv(1)), num2str(sf_uv(2)), num2str(sf_uv(3))];
    b_sync = '0';
    b_pitch = dec2bin(p_idx, 9);
    b_gain = dec2bin(g_idx-1, 10);
    b_rest = '00000000'; 
    
    bit_str = [b_sync, b_uv, '0', b_pitch, bits_lsf, b_gain, b_rest];
    % Pad/Truncate
    if length(bit_str) < 81, bit_str = [bit_str, repmat('0', 1, 81-length(bit_str))]; end
    all_bit_streams{s_idx} = bit_str(1:81);
    
    if mod(s_idx, 10) == 0, fprintf('.'); end
end
fprintf('\nDone.\n');
save('encoder_output.mat', 'all_bit_streams', 'signal_final');

% === MSVQ Helper ===
function [indices, q_vec] = run_msvq(target, cb_cells)
    resid = target;
    indices = zeros(1, 4);
    q_vec = zeros(size(target));
    
    for s = 1:4
        cb_stage = double(cb_cells{s});
        [~, idx] = min(sum((cb_stage - resid).^2, 1));
        indices(s) = idx;
        
        vec = cb_stage(:, idx);
        q_vec = q_vec + vec;
        resid = resid - vec;
    end
end

% === Pitch/Gain Placeholders (Please paste your functions here!) ===
function [p, v] = melp_pitch_integer(s, fs), p=50; v=0; end
function [p, v] = melp_pitch_fraction(s, p, fs), v=0; end
function [p, jit, v] = melp_pitch_final(s, fs, p, v), jit=0; end
% ... (其他函数)