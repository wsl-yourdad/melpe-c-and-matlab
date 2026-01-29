function [integer_pitch, corr_max] = melp_pitch_integer(input_frame, fs)
    % 整数基音周期提取 - 绝对安全版（防崩溃）
    
    % 默认返回值（清音情况）
    integer_pitch = 0;
    corr_max = 0;
    
    % 安全检查：信号太短直接返回
    if length(input_frame) < 160
        return;
    end
    
    % 1. 1KHz低通滤波
    Wn = 1000 / (fs/2);
    [b, a] = butter(4, Wn, 'low');
    filtered_frame = filter(b, a, input_frame);
    
    % 2. 计算自相关（40~160范围）
    min_lag = 40;
    max_lag = 160;
    [corr, lags] = xcorr(filtered_frame, max_lag, 'coeff');
    corr = corr(max_lag+1:end);
    lags = lags(max_lag+1:end);
    
    % 3. 取有效范围
    valid_mask = (lags >= min_lag) & (lags <= max_lag);
    valid_corr = corr(valid_mask);
    valid_lags = lags(valid_mask);
    
    % 4. 确保有有效值
    if isempty(valid_corr) || all(valid_corr == 0)
        return;  % 保持清音默认值
    end
    
    % 5. 找最大值
    [~, max_idx] = max(valid_corr);
    integer_pitch = valid_lags(max_idx);
    corr_max = valid_corr(max_idx);
    
    % 6. 清音检测（周期性弱）
    if corr_max < 0.3
        integer_pitch = 0;
        corr_max = 0;
    end
end