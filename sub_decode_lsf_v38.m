function lsf_mat = sub_decode_lsf_v38(p, prev_lsf, cb, res_cb_q15)
% =========================================================================
% V38 板块：与你“✅正常”的编码器完全匹配的反量化
% =========================================================================
    lsf_mat = zeros(10, 3);
    b = p.lsf_bits; 
    mv = double(cb.msvq_mean);
    
    if p.uv(3) == 0  % === Mode 1: UV ===
        % 1. F1/F2 绝对值 (编码器直接搜的 9bit)
        % 注意：编码器 pack 时减了1，这里 bin2dec 之后必须加1
        lsf_mat(:, 1) = sub_stabilize(double(cb.lsf_cb_9bit(:, bin2dec(b(1:9))+1)), 50);
        lsf_mat(:, 2) = sub_stabilize(double(cb.lsf_cb_9bit(:, bin2dec(b(10:18))+1)), 50);
        
        % 2. F3 MSVQ (8-6-5-5)
        m = [bin2dec(b(19:26))+1, bin2dec(b(27:32))+1, ...
             bin2dec(b(33:37))+1, bin2dec(b(38:42))+1];
        q_res = zeros(10, 1);
        for i = 1:4, q_res = q_res + double(cb.lsf_cb_mode_B{i}(:, m(i))); end
        lsf_mat(:, 3) = sub_stabilize(mv + q_res, 50);

    else  % === Mode 2: Voiced ===
        % 1. F3 Anchor (24 bits: 8-6-5-5)
        m = [bin2dec(b(1:8))+1, bin2dec(b(9:14))+1, ...
             bin2dec(b(15:19))+1, bin2dec(b(20:24))+1];
        q_res = zeros(10, 1);
        for i = 1:4, q_res = q_res + double(cb.lsf_cb_mode_B{i}(:, m(i))); end
        f3_q = sub_stabilize(mv + q_res, 50);
        lsf_mat(:, 3) = f3_q;
        
        % 2. 插值预测 (4 bits)
        int_idx = bin2dec(b(25:28)) + 1;
        w_vec = cb.inpCoef(int_idx, :); 
        w1 = w_vec(1:10)'; w2 = w_vec(11:20)';
        p1 = w1 .* prev_lsf + (1 - w1) .* f3_q;
        p2 = w2 .* prev_lsf + (1 - w2) .* f3_q;
        
        % 3. 残差 (14 bits: 8+6)
        r1 = bin2dec(b(29:36)) + 1;
        r2 = bin2dec(b(37:42)) + 1;
        % 这里的 res_cb_q15 必须是你在主脚本里除以 4 之后的
        res_20 = double(res_cb_q15{1}(:, r1)) + double(res_cb_q15{2}(:, r2));
        
        lsf_mat(:, 1) = sub_stabilize(p1 + res_20(1:10), 50);
        lsf_mat(:, 2) = sub_stabilize(p2 + res_20(11:20), 50);
    end
end

function lsf = sub_stabilize(lsf, min_d)
    lsf = sort(lsf(:));
    for i = 1:9
        if (lsf(i+1)-lsf(i)) < min_d
            lsf(i+1) = lsf(i) + min_d;
        end
    end
    lsf(lsf < 50) = 50; lsf(lsf > 32700) = 32700;
end