function disturbed_signal = add_interference(clean_signal, fs)
% 添加干扰 - 模拟论文图2.5
% 功能: 1. 添加50Hz单频干扰  2. 添加直流偏移
% 目的: 验证预处理滤波器能否干掉这些干扰
%
% 输入: clean_signal - 干净的语音信号
%       fs - 采样率
% 输出: disturbed_signal - 被"污染"的信号

% 1. 直流干扰（让信号整体抬高，不在0上下对称）
dc_offset = 0.1 * mean(abs(clean_signal));  % 抬高10%的平均幅度

% 2. 50Hz单频干扰（像交流电嗡嗡声）
t = (0:length(clean_signal)-1)' / fs;  % 时间轴
single_freq = 0.15 * max(abs(clean_signal)) * sin(2*pi*50*t);  % 50Hz正弦波

% 3. 加扰
disturbed_signal = clean_signal + dc_offset + single_freq;

% 4. 归一化（防止幅度过大溢出）
disturbed_signal = disturbed_signal / max(abs(disturbed_signal)) * 0.95;
end