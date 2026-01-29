function [pitches, bpvc_mat] = dequantize_pitch(p, cb)
    % 基于 uv_index 的 U/V 模式映射 (0->UUU, 1->VVU, 2->VUV, 4->UVV, 3/5/6/7->VVV)
    vo_map = [0 0 0; 1 1 0; 1 0 1; 1 1 1; 0 1 1; 1 1 1; 1 1 1; 1 1 1];
    uv_pattern = vo_map(p.uv_index + 1, :);
    
    pitches = zeros(1, 3);
    if p.uv_index == 0 % UUU
        pitches = [0 0 0];
    elseif any(p.uv_index == [1, 2, 4]) % UVV 等模式用 uvv_cb
        log_p = cb.pitch_vq_cb_uvv(:, p.pitch_idx + 1);
        pitches = 10.^log_p .* uv_pattern;
    else % VVV
        log_p = cb.pitch_vq_cb_vvv(:, p.pitch_idx + 1);
        pitches = 10.^log_p;
    end
    
    % 解码 BPVC (每帧 2 bits 映射到 5 频带)
    bpvc_mat = zeros(5, 3);
    for f = 1:3
        idx = bin2dec(p.bpvc_bits((f-1)*2+1 : f*2));
        % C代码标准：0->10000, 1->11000, 2->11100, 3->11111
        maps = [1 0 0 0 0; 1 1 0 0 0; 1 1 1 0 0; 1 1 1 1 1]';
        bpvc_mat(:, f) = maps(:, idx + 1) * uv_pattern(f);
    end
end