function [sig, st] = synthesis_filter(lsf, pitch, gc, gp, fsmag, st, df)
    FRAME = 180;
    % LSF -> LPC 弧度映射与转换
    a = lsf2poly((lsf(:)' / 32768) * pi);
    
    % 带宽扩展 gamma=0.994 (极点向内拉，防止发散爆鸣)
    a = a .* (0.994 .^ (0:10)); 
    
    % 稳定性防御：如果系数导致极点溢出，回退到上一帧
    if any(abs(roots(a)) >= 1.0)
        a = st.lpc_old; 
    else
        st.lpc_old = a;
    end
    
    % 激励生成与 Q15 能量补齐
    gt = linspace(gp, gc, FRAME);
    if pitch > 0
        T = round(pitch); if T < 20, T = 20; end
        pulse = zeros(1, FRAME);
        pulse(1:T:end) = sqrt(T) * 60; % 能量补偿系数，根据 Q15 调整
        [pulse, st.disp_mem] = filter(df, 1, pulse, st.disp_mem);
        exc = 0.8 * pulse + 0.2 * randn(1, FRAME) * 10;
    else
        exc = randn(1, FRAME) * 25; % 清音底噪
    end
    
    % 全反馈滤波 y(n) = x(n) - sum(a*y)
    [sig, st.syn_mem] = filter(1, a, exc .* gt, st.syn_mem);
    
    % 强制压限，防止产生 0dB 爆音
    sig(sig > 32767) = 32767; sig(sig < -32768) = -32768;
end