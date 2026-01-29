function [pitch_idx, pitch_quant] = melp_pitch_quantize(pitch, voicing)
% 功能: 将浮点基音周期量化为 7-bit 索引 (根据论文 2.2.4 节)
% 输入: pitch - 最终基音周期 (P3), voicing - 清浊音标志 (0 或 1)
% 输出: pitch_idx - 7-bit 量化索引 (0-127), pitch_quant - 量化后的基音值

    if voicing == 0 || pitch <= 0
        pitch_idx = 0; % 论文规则：清音帧量化为全零码字
        pitch_quant = 0;
        return;
    end

    % 1. 定义基音搜索范围 (采样点)
    P_MIN = 20; P_MAX = 160;
    
    % 2. 限制输入范围并取对数
    p_clamped = max(min(pitch, P_MAX), P_MIN);
    log_p = log10(p_clamped);
    
    % 3. 在对数域进行 99 电平均匀量化 (对应 20 到 160 采样点)
    log_min = log10(P_MIN);
    log_max = log10(P_MAX);
    
    % 计算量化索引 (映射到 1-99 范围)
    pitch_idx = round(1 + 98 * (log_p - log_min) / (log_max - log_min));
    
    % 4. 计算量化后的基音值 (反量化验证)
    pitch_quant = 10^(log_min + (pitch_idx - 1) * (log_max - log_min) / 98);
end                  