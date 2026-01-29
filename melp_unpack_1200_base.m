function [params_batch, next_prev_params] = melp_unpack_1200(bit_string, prev_params, codebooks)
% MELP_UNPACK_1200 (Final Strict Sync)
% 严格对照 pack_bits_c_style.m 逆向
% 修正了 LSF Mode 1/2 的分支读取逻辑

    %% 0. [码本字段映射]
    if ~isfield(codebooks, 'lsp_vq')
        if isfield(codebooks, 'lsf_cb_mode_C'), codebooks.lsp_vq = codebooks.lsf_cb_mode_C;
        elseif isfield(codebooks, 'lsf_cb'), codebooks.lsp_vq = codebooks.lsf_cb; end
    end
    if ~isfield(codebooks, 'gain') && isfield(codebooks, 'gain_cb'), codebooks.gain = codebooks.gain_cb; end
    if ~isfield(codebooks, 'fsvq_cb') && isfield(codebooks, 'fsvq_vq'), codebooks.fsvq_cb = codebooks.fsvq_vq; end
    if ~isfield(codebooks, 'pitch_vvv') && isfield(codebooks, 'pitch_vq_cb_vvv'), codebooks.pitch_vvv = codebooks.pitch_vq_cb_vvv; end

    %% 1. [比特流读取工具]
    if ischar(bit_string), bits = bit_string - '0'; else, bits = bit_string; end
    bit_ptr = 1;
    function val = read_bits(n_bits)
        if bit_ptr + n_bits - 1 > length(bits), val = 0; return; end
        chunk = bits(bit_ptr : bit_ptr + n_bits - 1);
        val = 0;
        for k = 1:n_bits, val = val * 2 + chunk(k); end
        bit_ptr = bit_ptr + n_bits;
    end

    %% 2. [解包流程]
    
    % (1) Sync (1 bit)
    sync_bit = read_bits(1);
    
    % (2) UV (3 bits)
    uv_raw = zeros(1,3);
    uv_raw(1) = read_bits(1);
    uv_raw(2) = read_bits(1);
    uv_raw(3) = read_bits(1);
    
    % (3) Pad (1 bit)
    pad_bit = read_bits(1);
    
    % (4) Pitch (9 bits) - 核心修正: Pitch 在 LSF 之前
    pitch_idx = read_bits(9);
    
    % (5) LSF (42 bits)
    % 判决: UV(3)=0 -> Mode 1 (UV); UV(3)=1 -> Mode 2 (Voiced)
    is_mode1 = (uv_raw(3) == 0); 
    
    lsf_f3 = prev_params.lsf; % 默认
    
    if is_mode1
        % === Mode 1 (UV): F1(9) + F2(9) + F3_MSVQ(24) ===
        idx_f1 = read_bits(9);
        idx_f2 = read_bits(9);
        
        % MSVQ (8+6+5+5 = 24)
        idx_f3_s1 = read_bits(8);
        idx_f3_s2 = read_bits(6);
        idx_f3_s3 = read_bits(5);
        idx_f3_s4 = read_bits(5);
        
        lsf_f3 = decode_msvq(codebooks, idx_f3_s1, idx_f3_s2, idx_f3_s3, idx_f3_s4);
        
        % 恢复 F1, F2 (使用 9bit 码本)
        if isfield(codebooks, 'lsf_cb_9bit')
            % 编码器打包做了 -1, 这里 idx 对应 MATLAB 索引要 +1
            lsf_f1 = codebooks.lsf_cb_9bit(:, min(idx_f1+1, size(codebooks.lsf_cb_9bit,2)));
            lsf_f2 = codebooks.lsf_cb_9bit(:, min(idx_f2+1, size(codebooks.lsf_cb_9bit,2)));
            % Q15 -> Rad 转换
            if max(abs(lsf_f1)) > 10, lsf_f1 = (lsf_f1/32768)*pi; end
            if max(abs(lsf_f2)) > 10, lsf_f2 = (lsf_f2/32768)*pi; end
        else
            lsf_f1 = lsf_f3; lsf_f2 = lsf_f3; % 降级
        end
        
    else
        % === Mode 2 (V): F3_MSVQ(24) + Int(4) + Res1(8) + Res2(6) ===
        idx_f3_s1 = read_bits(8);
        idx_f3_s2 = read_bits(6);
        idx_f3_s3 = read_bits(5);
        idx_f3_s4 = read_bits(5);
        
        lsf_f3 = decode_msvq(codebooks, idx_f3_s1, idx_f3_s2, idx_f3_s3, idx_f3_s4);
        
        int_idx = read_bits(4);
        res1 = read_bits(8); 
        res2 = read_bits(6);
        
        % 插值恢复 F1, F2
        w = int_idx / 15.0;
        lsf_f1 = (1-w)*prev_params.lsf + w*lsf_f3;
        lsf_f2 = (1-w)*prev_params.lsf + w*lsf_f3;
    end
    
    % (6) Gain (10 bits)
    gain_comb = read_bits(10);
    g_idx1 = floor(gain_comb / 32); 
    g_idx2 = mod(gain_comb, 32);
    
    % (7) BPVC (6 bits)
    bpvc_idx = read_bits(6);
    
    % (8) FSVQ (8 bits)
    fourier_idx = read_bits(8);
    
    % (9) Jitter (1 bit)
    jitter_idx = read_bits(1);
    
    %% 3. [参数还原]
    % Pitch
    min_p = 1.30103; max_p = 2.20412;
    step = (max_p - min_p) / 512;
    p_val = 10^(min_p + pitch_idx * step);
    if pitch_idx == 0, p_val = prev_params.pitch; end
    pitch_vec = [p_val; p_val; p_val];
    
    % Gain (自适应检查码本大小)
    if isfield(codebooks, 'gain')
        % 如果码本很大(>32)，说明是单级 VQ；如果是32，说明是 Split VQ
        if length(codebooks.gain) > 32
             % 可能是 10-bit VQ? 但1200bps通常不用这么大
             % 按照 5+5 split 处理最稳妥
             g1 = codebooks.gain(min(g_idx1+1, length(codebooks.gain)));
             g2 = codebooks.gain(min(g_idx2+1, length(codebooks.gain)));
        else
             g1 = codebooks.gain(min(g_idx1+1, length(codebooks.gain)));
             g2 = codebooks.gain(min(g_idx2+1, length(codebooks.gain)));
        end
        % Log -> Linear
        if g1 < 10, g1 = 10^g1; end
        if g2 < 10, g2 = 10^g2; end
        g1 = g1 * 200; g2 = g2 * 200;
    else
        g1=1000; g2=1000;
    end
    
    % BPVC
    bpvc_vec = zeros(5,1);
    for k=1:5, if bitand(bpvc_idx, 2^(k-1)), bpvc_vec(k)=1; end; end
    if pitch_idx > 0, bpvc_vec(1)=1; end
    
    % FSVQ
    fsvq_mags = ones(10,1);
    if isfield(codebooks, 'fsvq_cb')
        fsvq_mags = codebooks.fsvq_cb(:, min(fourier_idx+1, size(codebooks.fsvq_cb,2)));
        if max(fsvq_mags) > 5, fsvq_mags = fsvq_mags / max(fsvq_mags); end
    end
    
    %% 4. [输出组装]
    params_batch = repmat(prev_params, 3, 1);
    g_list = [g1, (g1+g2)/2, g2];
    
    for k=1:3
        params_batch(k).lsf = lpc_clamp(eval(['lsf_f',num2str(k)]), 0.005);
        params_batch(k).pitch = pitch_vec(k);
        params_batch(k).uv_flag = uv_raw(k);
        params_batch(k).gain = [g_list(k); g_list(k)];
        params_batch(k).jitter = jitter_idx;
        params_batch(k).bpvc = bpvc_vec;
        params_batch(k).fsvq = fsvq_mags;
    end
    next_prev_params = params_batch(3);
end

% 辅助: MSVQ求和
function lsf = decode_msvq(cb, i1, i2, i3, i4)
    % 索引保护
    i1 = min(i1+1, size(cb.lsp_vq{1},2));
    i2 = min(i2+1, size(cb.lsp_vq{2},2));
    i3 = min(i3+1, size(cb.lsp_vq{3},2));
    i4 = min(i4+1, size(cb.lsp_vq{4},2));
    
    lsf = cb.lsp_vq{1}(:, i1) + cb.lsp_vq{2}(:, i2) + ...
          cb.lsp_vq{3}(:, i3) + cb.lsp_vq{4}(:, i4);
          
    if max(abs(lsf)) > 10, lsf = (lsf / 32768.0) * pi; end
end

function lsf = lpc_clamp(lsf, delta)
    lsf = sort(lsf);
    lsf(lsf<delta)=delta; lsf(lsf>pi-delta)=pi-delta;
end