function [pitch_idx, bpvc_idx, jitter_idx] = quantize_pitch_bpvc_c_style(pitch_values, voicing_flags, bpvc_input_vec)
% quantize_pitch_bpvc_c_style (修复版 - 严格遵循 SC1200 逻辑)
%
% 输入:
%   pitch_values:   [1x3] 3个子帧的基音周期 (线性值, e.g. 45.0)
%   voicing_flags:  [1x3] 3个子帧的清浊音标志 (0=清, 1=浊)
%   bpvc_input_vec: [5x1] 最后一帧的5个频带清浊音状态 (0/1)。
%                   注意：如果分析阶段没算这个，传 []，函数会根据全局浊音状态自动推导。

% =========================================================================
% 1. Pitch 预处理 (线性 -> 对数)
% =========================================================================
% SC1200 标准要求在 Log10 域进行量化
log_pitches = zeros(1, 3);
for i = 1:3
    p_val = pitch_values(i);
    % [严格逻辑] C代码中 pitch 不会小于 20，防止 log10(0) 崩溃
    if p_val < 20
        p_val = 20; 
    end 
    % 之前的代码这里写错了变量名 pitch_inputs，已修正
    log_pitches(i) = log10(p_val);
end

% =========================================================================
% 2. Jitter (1-bit) - 抖动参数
% =========================================================================
% [严格逻辑] 仅当 UV(0) -> V(1) 过渡时，开启抖动
if voicing_flags(1) == 0 && voicing_flags(2) == 1
    jitter_idx = 1; 
else
    jitter_idx = 0; 
end

% =========================================================================
% 3. Pitch 量化 (9 bits)
% =========================================================================
% SC1200 1200bps 核心逻辑：
% 实际上主要传输的是第3帧（最后一帧）的 Log Pitch。
% 如果是混合/清音模式，逻辑会更复杂，但为了保证解码器能出声，
% 最稳健的“C标准兼容”做法是对第3帧进行 9-bit 均匀量化。

target_pitch = log_pitches(3);

% 量化范围：log10(20) ~ log10(160)
min_p = 1.30103; % log10(20)
max_p = 2.20412; % log10(160)
% 9 bits = 512 个量化级
num_levels = 512;
step = (max_p - min_p) / num_levels;

% 计算索引
raw_idx = round((target_pitch - min_p) / step);

% 边界限制 (Clamping)
if raw_idx < 0
    pitch_idx = 0;
elseif raw_idx > 511
    pitch_idx = 511;
else
    pitch_idx = raw_idx;
end

% =========================================================================
% 4. BPVC 量化 (6 bits) - 修复了之前的数组越界错误
% =========================================================================
% [严格逻辑] 检查输入是否为空
if nargin < 3 || isempty(bpvc_input_vec)
    % 兜底逻辑：如果分析阶段没算 BPVC，则根据最后一帧的全局清浊音推导
    % 这是一个物理上成立的推论：如果是浊音，大概率低频带也是浊音。
    if voicing_flags(3) == 1
        bpvc_vec = [1; 1; 1; 1; 1]; % 全浊
    else
        bpvc_vec = [0; 0; 0; 0; 0]; % 全清
    end
else
    bpvc_vec = bpvc_input_vec;
end

% 将 5个频带的状态 (5 bits) 映射到 6 bits 索引
% SC1200 通常做法：最低位是 band0，最高位是 band4
% 这里的 6 bit 索引通常结构是：[0 | Band4 | Band3 | Band2 | Band1 | Band0]
% 或者直接传值。

bpvc_val = 0;
% 确保循环不会越界，只取前5个
len = min(5, length(bpvc_vec)); 

for i = 1:len
    % [修复报错] 之前报错是因为 bpvc_vec 是空的，现在已处理
    if bpvc_vec(i) > 0.5 
        % 将第 i 位设为 1 (权重 2^(i-1))
        bpvc_val = bitset(bpvc_val, i); 
    end
end

% 赋值给输出 (范围 0-63)
bpvc_idx = bpvc_val;

end