function [fractional_pitch, Vp] = melp_pitch_fraction(input_frame, integer_pitch, fs)
    % 分数基音周期提取 - 安全版
    
    % 默认返回值
    fractional_pitch = integer_pitch;  % 失败时用整数周期
    Vp = 0;
    
    % 安全检查
    if integer_pitch == 0 || length(input_frame) < 160
        return;
    end
    
    % 1. 500Hz子带滤波
    Wn = 500 / (fs/2);
    [b, a] = butter(4, Wn, 'low');
    subband_frame = filter(b, a, input_frame);
    
    % 2. 自相关计算
    search_range = 20:160;
    [corr, ~] = xcorr(subband_frame, max(search_range), 'coeff');
    corr = corr(max(search_range)+1:end);
    
    % 3. 插值
    T = integer_pitch;
    if T > 20 && T < 160 && T+1 <= length(corr)
        r_tm1 = corr(T-1);
        r_t = corr(T);
        r_tp1 = corr(T+1);
        
        % 抛物线插值
        denom = r_tm1 + r_tp1 - 2*r_t;
        if denom ~= 0
            delta = (r_tm1 - r_tp1) / (2*denom);
            delta = max(min(delta, 1), -1);
            fractional_pitch = T + delta;
            
            % Vp计算
            Vp = r_t - 0.25 * (r_tm1 - r_tp1) * delta;
            Vp = max(min(Vp, 1), 0);
        else
            fractional_pitch = T;
            Vp = r_t;
        end
    else
        fractional_pitch = T;
        Vp = corr(T);
    end
end