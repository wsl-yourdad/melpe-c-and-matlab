function [frame_sig, st] = synthesis_engine(lsf_s, lsf_e, p_s, p_e, bpvc, g_s, g_e, fsmag, st, df, bpf_n, bpf_d)
% synthesis_engine: 12个输入参数
% =========================================================================
    FRAME_LEN = 180;
    NUM_SUB = 4;
    SUB_LEN = FRAME_LEN / NUM_SUB;
    frame_sig = zeros(1, FRAME_LEN);
    
    % 1. 动态重塑 BPF 码本 (45x1 -> 5x9)
    bpf_num_all = reshape(bpf_n, 9, 5)';
    bpf_den_all = reshape(bpf_d, 9, 5)';

    % 2. 准备线性轨迹
    g_traj = linspace(g_s, g_e, FRAME_LEN);
    p_traj = linspace(p_s, p_e, FRAME_LEN);

    for i = 1:NUM_SUB
        idx_range = (i-1)*SUB_LEN + 1 : i*SUB_LEN;
        alpha = (i - 0.5) / NUM_SUB;
        
        % A. LSF 插值与 A(z) 稳定性
        lsf_sub = (1-alpha)*lsf_s + alpha*lsf_e;
        a = lsf2poly((lsf_sub' / 32768) * pi);
        a = a .* (0.994 .^ (0:10)); 
        if any(abs(roots(a)) >= 1.0), a = st.lpc_old; else, st.lpc_old = a; end

        % B. 混合激励源构建
        curr_p = p_traj((i-1)*SUB_LEN + 1);
        noise = randn(1, SUB_LEN);
        pulse = zeros(1, SUB_LEN);
        
        if curr_p > 20
            T = round(curr_p);
            pulse(mod(idx_range, T) == 0) = sqrt(T) * 40; 
            % 傅里叶增强
            w0 = 2 * pi / curr_p;
            enh = zeros(1, SUB_LEN);
            for k = 1:10
                enh = enh + (fsmag(k)/8192) * cos(w0 * (idx_range-1) + st.f_phase(k));
                st.f_phase(k) = mod(st.f_phase(k) + w0 * SUB_LEN, 2*pi);
            end
            pulse = pulse + enh * 10;
            [pulse, st.disp_mem] = filter(df, 1, pulse, st.disp_mem);
        end

        % C. 5-Band 混合 (真正用到 45x1 码本的地方)
        mixed_exc = zeros(1, SUB_LEN);
        for b = 1:5
            if bpvc(b) == 1, src = pulse; else, src = noise * 15; end
            [band_out, st.bpf_mem(b,:)] = filter(bpf_num_all(b,:), bpf_den_all(b,:), src, st.bpf_mem(b,:));
            mixed_exc = mixed_exc + band_out;
        end
        
        % D. LPC 综合
        [sub_sig, st.syn_mem] = filter(1, a, mixed_exc, st.syn_mem);
        lin_g = 10^(g_traj(idx_range(1))/20);
        frame_sig(idx_range) = sub_sig * lin_g;
    end
end