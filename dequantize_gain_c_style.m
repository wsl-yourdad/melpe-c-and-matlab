function gains_dequant = dequantize_gain_c_style(gain_idx, gain_codebook)
% dequantize_gain_c_style: 根据10-bit索引，从增益码本中查表
%
% 输入:
%   gain_idx:       10-bit 增益索引 (1-1024)
%   gain_codebook:  1024x6 的增益码本
%
% 输出:
%   gains_dequant:  [2x3] 反量化后的6个增益值
%

% C语言的10-bit VQ是直接查表
if gain_idx > 0 && gain_idx <= size(gain_codebook, 1)
    gains_6D = gain_codebook(gain_idx, :);
else
    % 如果索引无效，返回一个默认的低能量增益
    gains_6D = ones(1, 6) * 20; % 20dB
end

% 将6维向量重塑为2x3矩阵
gains_dequant = reshape(gains_6D, 2, 3);

end
