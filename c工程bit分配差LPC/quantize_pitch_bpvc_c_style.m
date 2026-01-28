function [pitch_idx, bpvc_idx, jitter_idx] = quantize_pitch_bpvc_c_style(pitch_values, voicing_flags, bpvc_input_vec)
% quantize_pitch_bpvc_c_style (防弹版)
% 修复内容：强制初始化 BPVC 向量，防止空输入导致索引越界

% =========================================================================
% 1. Pitch 预处理 (线性 -> 对数)
% =========================================================================
log_pitches = zeros(1, 3);
for i = 1:3
    p_val = pitch_values(i);
    if p_val < 20, p_val = 20; end 
    log_pitches(i) = log10(p_val);
end

% =========================================================================
% 2. Jitter (1-bit)
% =========================================================================
if voicing_flags(1) == 0 && voicing_flags(2) == 1
    jitter_idx = 1; 
else
    jitter_idx = 0; 
end

% =========================================================================
% 3. Pitch 量化 (9 bits)
% =========================================================================
target_pitch = log_pitches(3);
min_p = 1.30103; % log10(20)
max_p = 2.20412; % log10(160)
step = (max_p - min_p) / 512;
raw_idx = round((target_pitch - min_p) / step);

if raw_idx < 0
    pitch_idx = 0;
elseif raw_idx > 511
    pitch_idx = 511;
else
    pitch_idx = raw_idx;
end

% =========================================================================
% 4. BPVC 量化 (6 bits) - 【高保真修正版】
% =========================================================================
% 输入检查：严禁为空！必须由分析器提供真实的 5维 BPVC 向量
if nargin < 3 || isempty(bpvc_input_vec)
    % 如果是清音帧，BPVC 确实可以是全0，但必须显式传入
    if voicing_flags(3) == 0
        bpvc_vec = zeros(5,1);
    else
        % 既然是浊音，就必须有真实的带通分析数据，不能瞎猜！
        error('❌ 严重逻辑错误：Frame 3 为浊音，但 BPVC 输入向量为空！请检查 melp_ana 分析模块。');
    end
else
    % 直接使用分析模块算出来的真实值 (可能是 [1;1;0;0;0] 这种混合状态)
    bpvc_vec = bpvc_input_vec; 
end

% 接下来进行正常的 BPVC 编码 (通常是映射到 60 个合法状态之一或直接量化)
% 转换逻辑 (5 bits -> 整数)
bpvc_val = 0;
% 动态获取长度，防止越界 (即使 bpvc_vec 是 5x1，这里也稳健)
safe_len = min(5, length(bpvc_vec)); 

for i = 1:safe_len
    if bpvc_vec(i) > 0.5 
        bpvc_val = bitset(bpvc_val, i); 
    end
end

bpvc_idx = bpvc_val;

end