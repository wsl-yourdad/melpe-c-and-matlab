function [pitches, bpvc_mat] = dequantize_pitch_joint(p, cb)
    % p.uv_index: 3-bit 决定超帧的清浊音模式 (0-7)
    % p.pitch_idx: 9-bit 基音索引
    
    % vo_map 严格对照 C 代码 1200bps 模式映射
    % 0:UUU, 1:VVU, 2:VUV, 4:UVV, 3/5/6/7: 不同偏移下的 VVV
    vo_map = [0 0 0; 1 1 0; 1 0 1; 1 1 1; 0 1 1; 1 1 1; 1 1 1; 1 1 1];
    uv_pattern = vo_map(p.uv_index + 1, :);
    
    pitches = zeros(1, 3);
    
    if p.uv_index == 0 % 全清音 UUU
        pitches = [0 0 0];
    elseif any(p.uv_index == [1, 2, 4]) % 只有两帧是浊音
        % 使用 pitch_vq_cb_uvv 码本
        log_p = cb.pitch_vq_cb_uvv(:, p.pitch_idx + 1);
        pitches = (10.^log_p) .* uv_pattern'; % 只有 V 的帧有基音
    else % VVV 模式
        % 使用 pitch_vq_cb_vvv 码本
        log_p = cb.pitch_vq_cb_vvv(:, p.pitch_idx + 1);
        pitches = 10.^log_p;
    end
    
    % 解码 BPVC (Bandpass Voicing)
    % C代码: 每帧 2bits 映射到 5 个频带的清浊音标志
    bpvc_mat = zeros(5, 3);
    % 解包 bpvc_bits (6位 -> 每帧2位)
    for f = 1:3
        idx = bin2dec(p.bpvc_bits((f-1)*2+1 : f*2));
        % C代码标准映射表: 0->低频1, 1->2, 2->3, 3->全V
        maps = [1 0 0 0 0; 1 1 0 0 0; 1 1 1 0 0; 1 1 1 1 1]';
        bpvc_mat(:, f) = maps(:, idx + 1) * uv_pattern(f);
    end
end