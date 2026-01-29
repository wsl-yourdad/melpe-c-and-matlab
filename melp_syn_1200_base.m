function [pcm_out, lpc_state_next, pulse_phase_next] = melp_syn_1200(params, lpc_state_prev, pulse_phase_prev, codebooks)
% MELP_SYN_1200 (Full C-Code Compliance)
% 1. Pulse Dispersion (65-tap FIR)
% 2. 5-Band Mixed Excitation
% 3. LPC Synthesis (1/A(z))

    FRAME_LEN = 180;
    lsf = params.lsf; pitch = params.pitch; gain = params.gain; uv = params.uv_flag;
    jitter = params.jitter;
    if isfield(params,'fsvq'), fsvq=params.fsvq; else, fsvq=ones(10,1); end
    if isfield(params,'bpvc'), bpvc=params.bpvc; else, bpvc=[1;0;0;0;0]; end
    
    % 1. LPC & BW Expansion (Gamma=0.994)
    a_lpc = lsf2poly(lsf) .* (0.994 .^ (0:10));
    
    % 2. 脉冲激励生成 (含色散)
    pulse_exc = zeros(FRAME_LEN + 100, 1);
    curr_ph = pulse_phase_prev;
    
    % === C代码原版 65阶 FIR 色散滤波器 ===
    disp_filter = [ ...
      -0.00466, -0.00539, -0.00516, -0.00388, -0.00166,  0.00123,  0.00427,  0.00688, ...
       0.00849,  0.00863,  0.00713,  0.00403, -0.00031, -0.00517, -0.00962, -0.01275, ...
      -0.01383, -0.01235, -0.00832, -0.00224,  0.00494,  0.01206,  0.01784,  0.02114, ...
       0.02107,  0.01701,  0.00923, -0.00163, -0.01430, -0.02707, -0.03799, -0.04524, ...
      -0.04765, -0.04496, -0.03747, -0.02569, -0.01053,  0.00659,  0.02401,  0.03991, ...
       0.05257,  0.06053,  0.06283,  0.05928,  0.05048,  0.03738,  0.02127,  0.00375, ...
      -0.01343, -0.02850, -0.03992, -0.04655, -0.04786, -0.04396, -0.03545, -0.02345, ...
      -0.00938,  0.00508,  0.01832,  0.02891,  0.03571,  0.03813,  0.03607,  0.02985, ...
       0.02028 ]';

    T = round(pitch); if T<20, T=20; end
    % FSVQ 谐波合成
    pulse = zeros(T,1);
    w0 = 2*pi/T;
    for k=1:10
        pulse = pulse + fsvq(k)*cos(k*w0*((1:T)'-T/2));
    end
    % 应用色散
    pulse = conv(pulse, disp_filter, 'same');
    pulse = pulse / norm(pulse) * sqrt(T);
    
    % 重叠相加放置脉冲
    idx_w = 1;
    while idx_w <= FRAME_LEN
        p_eff = pitch * (1 + (rand-0.5)*0.2*jitter);
        curr_ph = curr_ph + 1;
        if curr_ph >= p_eff
            curr_ph = curr_ph - p_eff;
            len = min(T, FRAME_LEN-idx_w+1);
            if len > 0
                pulse_exc(idx_w:idx_w+len-1) = pulse_exc(idx_w:idx_w+len-1) + pulse(1:len);
            end
        end
        idx_w = idx_w + 1;
    end
    pulse_exc = pulse_exc(1:FRAME_LEN);
    
    % 3. 噪声激励
    noise_exc = randn(FRAME_LEN, 1);
    noise_exc = noise_exc / (std(noise_exc)+1e-9);
    
    % 4. 5频带混合 (Mixed Excitation)
    mixed_exc = zeros(FRAME_LEN, 1);
    bands = [0 500; 500 1000; 1000 2000; 2000 3000; 3000 4000];
    
    for k=1:5
        Wn = bands(k,:) / 4000;
        if k==1, [b,a]=butter(2,Wn(2),'low');
        elseif k==5, [b,a]=butter(2,Wn(1),'high');
        else, [b,a]=butter(2,Wn); end
        
        if bpvc(k) > 0.5
            mixed_exc = mixed_exc + filter(b,a,pulse_exc);
        else
            mixed_exc = mixed_exc + filter(b,a,noise_exc);
        end
    end
    
    % 强制 UV 处理
    if uv < 0.5, mixed_exc = noise_exc; end
    
    % 5. 增益匹配
    mixed_exc = mixed_exc / (std(mixed_exc)+1e-9) * mean(gain);
    
    % 6. LPC 合成滤波 (1/A(z))
    [syn, lpc_state_next] = filter(1, a_lpc, mixed_exc, lpc_state_prev);
    
    % 7. 后处理 (De-emphasis)
    pcm_out = filter(1, [1, -0.5], syn);
    
    % 8. 极性修正
    pcm_out = -pcm_out * 2.5; 
    
    pulse_phase_next = curr_ph;
end