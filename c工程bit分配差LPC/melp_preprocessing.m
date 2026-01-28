function filtered_signal = melp_preprocessing(input_signal, fs)
% MELP预处理滤波器 - 论文公式(2-3)实现
% 功能: 去除直流分量和60Hz以下低频噪声
% 滤波器: 四阶巴特沃斯，截止频率60Hz，阻带30dB
% 参考文献: 论文第13页
% 
% 输入: 
%   input_signal - 原始语音信号
%   fs - 采样率(Hz)
% 输出:
%   filtered_signal - 预处理后的信号

% 归一化截止频率（60Hz @ 8kHz）
Wn = 60 / (fs/2); 

% 四阶巴特沃斯滤波器系数（论文公式2-3的MATLAB实现）
% 注意: 论文给的是传递函数分母，butter函数会自动计算
[b, a] = butter(4, Wn, 'high');  % 'high'是高通，因为我们想滤掉60Hz以下

% 应用滤波器
filtered_signal = filter(b, a, input_signal);

end