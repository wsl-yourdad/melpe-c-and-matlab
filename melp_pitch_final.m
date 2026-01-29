function [final_pitch, Vp_final, Voicing] = melp_pitch_final(residual_frame, fs, pitch_candidate, Vp_subband)
    % 最终基音周期计算 - 增强型保护版
    final_pitch = 0; Vp_final = 0; Voicing = 0;
    
    % 1. 预处理：1KHz低通滤波
    Wn = 1000 / (fs/2); [b, a] = butter(6, Wn, 'low');
    residual_filtered = filter(b, a, residual_frame);
    
    % 2. 计算残差相关性 (p3_int 可能是整数)
    [p3_int, Vp_p3] = melp_pitch_integer(residual_filtered, fs);
    
    % 3. 维护长时平均 Pavg
    persistent pitch_history;
    if isempty(pitch_history) || length(pitch_history) < 3
        pitch_history = [60, 60, 60];
    end
    Pavg = mean(pitch_history);
    
    % 4. 逻辑判决：参考原始信号强度 Vp_subband
    % 💡 核心逻辑：如果原始信号相关性(Vp_subband)很高，即使残差(Vp_p3)很低，也要维持浊音
    if Vp_p3 >= 0.6
        Dth = (p3_int <= 100) * 0.75 + (p3_int > 100) * 0.5;
        [best_int, Vp_res] = melp_multipulse_detect(residual_filtered, fs, p3_int, Vp_p3, Dth);
        Vp_final = Vp_res;
    else
        % 弱浊音/清音保护分支
        if Vp_p3 < 0.55 && Vp_subband < 0.6 % 真正判定为清音
            best_int = round(Pavg);
            Vp_final = Vp_p3;
        else
            % 如果原始信号强(Vp_subband > 0.6)，即便残差弱，也认为是浊音
            Dth = (p3_int <= 100) * 0.9 + (p3_int > 100) * 0.7;
            [best_int, Vp_res] = melp_multipulse_detect(residual_filtered, fs, p3_int, max(Vp_p3, 0.55), Dth);
            Vp_final = max(Vp_res, Vp_subband * 0.85); % 借用原始信号的相关性
        end
    end
    
    % 5. 分数插值与最终输出
    if best_int > 0
        [refined_pitch, ~] = melp_pitch_fraction(residual_filtered, round(best_int), fs);
    else
        refined_pitch = 0;
    end
    
    % 设置更稳健的门限
    if Vp_final >= 0.55 && refined_pitch > 0 %改为>=0.45后第42帧的最终基音周期将变为非0值
        %同时，整个5s样本内的浊音帧的个数也会增多
        final_pitch = refined_pitch;
        Voicing = 1;
        pitch_history = [pitch_history(2:end), final_pitch];
    else
        final_pitch = 0; Voicing = 0;
    end
end