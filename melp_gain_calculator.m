function [G1, G2] = melp_gain_calculator(speech_frame, residual_frame, pitch, Vp_final)
    % [G1, G2] = melp_gain_calculator(...)
    % 修正版：强制输出 dB 值以匹配 C 码本 (10.0 ~ 74.0)
    
    % --- 步骤 1: 确定窗口大小 L (逻辑保持不变) ---
    if Vp_final < 0.6
        L = 120;
    else
        P2 = round(pitch);
        if P2 == 0
            L = 120;
        else
            n = 1;
            L = n * P2;
            while L <= 120
                n = n + 1;
                L = n * P2;
            end
            if L > 320
                L = L / 2;
            end
        end
    end
    L = round(L); 
    
    % --- 步骤 2: 计算 G2 (语音增益) ---
    ref_point_idx = length(speech_frame); 
    start_G2 = max(1, ref_point_idx - L + 1);
    end_G2 = ref_point_idx;
    window_G2 = speech_frame(start_G2:end_G2);
    
    % 计算线性 RMS
    g2_linear = melp_calculate_rms(window_G2);
    
    % --- 步骤 3: 计算 G1 (残差增益) ---
    ref_G1 = ref_point_idx - 9;
    start_G1 = max(1, ref_G1 - L + 1);
    end_G1 = ref_G1;
    
    if end_G1 < 1
        window_G1 = []; 
        g1_linear = 0.01; % 防止空窗口
    else
        window_G1 = residual_frame(start_G1:end_G1);
        g1_linear = melp_calculate_rms(window_G1);
    end
    
    % =====================================================================
    % [核心修正] 线性幅度 -> 分贝 (dB)
    % =====================================================================
    % 码本范围是 10~74。
    % 假设输入信号是 Q15 (±32767)，那么 RMS 最大约 20000。
    % 20*log10(20000) ≈ 86 dB。
    % 20*log10(10) ≈ 20 dB。
    % 这与码本完美匹配！
    
    % 增加 0.01 防止 log10(0)
    G1 = 20 * log10(g1_linear + 0.01);
    G2 = 20 * log10(g2_linear + 0.01);
    
    % 钳位保护 (可选，防止极小噪声导致 -Inf)
    if G1 < 0, G1 = 0; end
    if G2 < 0, G2 = 0; end
end

% --- 内部辅助函数: 计算 RMS ---
% (我直接把它写在这里，防止你原来的 melp_calculate_gain 内部有其他缩放)
function rms_val = melp_calculate_rms(vec)
    if isempty(vec)
        rms_val = 0.01;
        return;
    end
    energy = sum(double(vec).^2); % 转double防止溢出
    rms_val = sqrt(energy / length(vec));
end