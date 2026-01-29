function gains_dequant = dequantize_gain(gain_idx, gain_cb_matrix)
% dequantize_gain: 直接对增益码本矩阵进行行查表
% gain_idx: 10-bit 索引 (0-1023)
% gain_cb_matrix: 已经在主脚本传入的 CB.gain 矩阵

    idx = gain_idx + 1; % 转为 MATLAB 1-based 索引
    
    % 直接对传入的矩阵 gain_cb_matrix 进行尺寸检查
    if idx > 0 && idx <= size(gain_cb_matrix, 1)
        gains_6D = gain_cb_matrix(idx, :);
    else
        % 索引越界保护，给个默认 dB 值
        gains_6D = ones(1, 6) * 10; 
    end
    
    % 重塑为 2x3 结构 (每帧 2 个子帧增益，共 3 帧)
    gains_dequant = reshape(gains_6D, 2, 3);
end