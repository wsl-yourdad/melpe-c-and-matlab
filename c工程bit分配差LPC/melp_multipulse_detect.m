function [detected_pitch, Vp_detected] = melp_multipulse_detect(residual_frame, fs, pitch_candidate, Vp_input, Dth)
    % 1. 安全检查：如果候选值为 0，直接返回
    if pitch_candidate <= 0
        detected_pitch = 0; Vp_detected = 0; return;
    end

    % 2. 约束搜索范围
    min_search = max(round(pitch_candidate - 5), 20);
    max_search = min(round(pitch_candidate + 5), 160);
    
    % 3. 局部搜索 (修正索引逻辑)
    [corr, lags] = xcorr(residual_frame, 160, 'coeff');
    center_idx = 161;
    corr_section = corr(center_idx + (min_search:max_search));
    lags_section = min_search : max_search;
    
    [Vp_int, max_idx] = max(corr_section);
    p_int = lags_section(max_idx);
    
    % 4. 判决
    if Vp_int > Vp_input * Dth
        detected_pitch = p_int; Vp_detected = Vp_int;
    else
        detected_pitch = pitch_candidate; Vp_detected = Vp_input;
    end
end