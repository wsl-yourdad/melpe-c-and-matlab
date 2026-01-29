function lsf_mat = dequantize_lsf(p, prev_lsf, cb, res_raw)
    % 严格按照 SC1200 low_rate_chn_read 逻辑
    lsf_mat = zeros(10, 3);
    b = p.lsf_bits; 
    mv = double(cb.msvq_mean);
    
    % --- 确定三帧 LSF (F3 是锚点) ---
    if p.uv_index == 0 % UUU
        lsf_mat(:,1) = lsf_stabilize(double(cb.lsf_cb_9bit(:, bin2dec(b(1:9))+1)));
        lsf_mat(:,2) = lsf_stabilize(double(cb.lsf_cb_9bit(:, bin2dec(b(10:18))+1)));
        m = [bin2dec(b(19:26)), bin2dec(b(27:32)), bin2dec(b(33:37)), bin2dec(b(38:42))] + 1;
        q3 = mv; for i=1:4, q3 = q3 + double(cb.lsf_cb_mode_B{i}(:,m(i))); end
        lsf_mat(:,3) = lsf_stabilize(q3);
    else
        % F3 Anchor (24 bits)
        m = [bin2dec(b(1:8)), bin2dec(b(9:14)), bin2dec(b(15:19)), bin2dec(b(20:24))] + 1;
        q3_res = zeros(10, 1);
        for i=1:4, q3_res = q3_res + double(cb.lsf_cb_mode_B{i}(:, m(i))); end
        lsf_f3 = lsf_stabilize(mv + q3_res);
        lsf_mat(:, 3) = lsf_f3;

        % F1/F2 Prediction (Q14 Weights)
        int_idx = bin2dec(b(25:28)) + 1;
        w = cb.inpCoef(int_idx, :); 
        p1 = (w(1:10)' .* prev_lsf + (16384-w(1:10)') .* lsf_f3) / 16384;
        p2 = (w(11:20)' .* prev_lsf + (16384-w(11:20)') .* lsf_f3) / 16384;
        
        % F1/F2 Residuals (20-Dim)
        r = [bin2dec(b(29:36)), bin2dec(b(37:42))] + 1;
        res_20 = double(res_raw(:, r(1))) + double(res_raw(:, 256 + r(2)));
        
        % Final Addition (Residual / 4 as per C code)
        lsf_mat(:,1) = lsf_stabilize(p1 + res_20(1:10)/4.0);
        lsf_mat(:,2) = lsf_stabilize(p2 + res_20(11:20)/4.0);
    end
end

function lsf = lsf_stabilize(lsf)
    lsf = sort(lsf(:));
    min_dist = 60;
    for i = 1:9
        if (lsf(i+1)-lsf(i)) < min_dist, lsf(i+1) = lsf(i) + min_dist; end
    end
    lsf(lsf<100)=100; lsf(lsf>32600)=32600;
end