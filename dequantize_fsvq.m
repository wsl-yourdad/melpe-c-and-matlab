function fsmag = dequantize_fsvq(fs_idx, cb)
    % fourier_idx: 8-bit 索引 (0-255)
    % cb.fsvq_cb: [10 x 256] 傅里叶幅值码本
    
    % MATLAB 索引从 1 开始
    idx = fs_idx + 1;
    
    if idx > 0 && idx <= size(cb.fsvq_cb, 2)
        fsmag = double(cb.fsvq_cb(:, idx));
    else
        % 默认平滑幅值 (C代码标准中的 ONE_Q13)
        fsmag = ones(10, 1) * 8192; 
    end
end