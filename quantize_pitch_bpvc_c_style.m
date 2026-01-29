function [pitch_idx_out, bpvc_idx, jitter_idx, uv_encoded_idx] = quantize_pitch_bpvc_c_style(pitch_values, voicing_flags, bpvc_input_vec, codebooks)
% =========================================================================
% Pitch VQ (Vector Quantization) 1200bps C-Style - FIXED
% 修正: 增加 1200bps 特有的 Pitch/UV 联合编码逻辑
%       (VVV模式下将 11bit Pitch 拆分为 9bit Pitch + 3bit Encoded UV)
% =========================================================================

    % 1. 预处理: Log10 域 + 钳位
    log_pitches = zeros(3, 1);
    for i = 1:3
        p = pitch_values(i);
        if p < 20, p = 20; end
        if p > 160, p = 160; end
        log_pitches(i) = log10(p);
    end
    
    UV_PITCH_LOG = log10(50.0);

    % 2. 准备 VQ 目标向量与权重
    target = log_pitches;
    weights = double(voicing_flags(:)); % 1=Voiced
    
    % 处理清音帧的目标值
    for i = 1:3
        if weights(i) == 0 % UV
            target(i) = UV_PITCH_LOG;
        end
    end
    
    cnt_voiced = sum(weights);
    
    % 3. 矢量量化逻辑 (Match pitch_vq)
    % ---------------------------------------------------------------------
    if cnt_voiced <= 1
        % 全清音或微弱浊音 -> 使用 UVV 码本
        cb_use = codebooks.pitch_vq_cb_uvv; % 512 entries
    else
        if cnt_voiced == 3
            cb_use = codebooks.pitch_vq_cb_vvv; % 2048 entries
        else
            cb_use = codebooks.pitch_vq_cb_uvv; % 512 entries
        end
    end
    
    % 执行 VQ 搜索
    best_err = inf;
    best_idx = 1;
    [dim, cb_size] = size(cb_use);
    
    for i = 1:cb_size
        diff = target - cb_use(:, i);
        err = sum(weights .* (diff.^2));
        if err < best_err
            best_err = err;
            best_idx = i;
        end
    end
    
    raw_pitch_idx = best_idx - 1; % 0-based index

    % =====================================================================
    % [核心修正] 1200bps 编码逻辑: Pitch/UV 联合映射
    % 参考 C代码 low_rate_chn_write 中的逻辑
    % =====================================================================
    
    % VVV Index Mapping Table (from qnt12_cb.c)
    % Maps high bits 0,1,2,3 -> UV Codes 3,5,6,7
    vvv_map = [3, 5, 6, 7]; 
    
    if cnt_voiced == 3
        % === VVV Mode (All Voiced) ===
        % Input Pitch Index: 0 ~ 2047 (11 bits)
        % Split into: High 2 bits (-> UV field) + Low 9 bits (-> Pitch field)
        
        idx_high = floor(raw_pitch_idx / 512); % 取高2位 (0-3)
        pitch_idx_out = mod(raw_pitch_idx, 512); % 取低9位
        
        % Map high bits to UV Code
        uv_encoded_idx = vvv_map(idx_high + 1);
        
    elseif cnt_voiced == 2
        % === Mixed Mode (UVV, VUV, VVU) ===
        % Pitch Index: 0 ~ 511 (9 bits) - No change
        pitch_idx_out = raw_pitch_idx;
        
        % Map Flags to UV Code (C-Style Logic)
        % 1=Voiced in MATLAB.
        if voicing_flags(1) == 0     % U V V
            uv_encoded_idx = 4;
        elseif voicing_flags(2) == 0 % V U V
            uv_encoded_idx = 2;
        else                         % V V U (implies flags(3)==0)
            uv_encoded_idx = 1;
        end
        
    else
        % === Mostly UV (UUU, UUV, UVU, VUU) ===
        % Pitch Index: 0 ~ 511
        pitch_idx_out = raw_pitch_idx;
        uv_encoded_idx = 0; 
        
        % Note: C code might remap pitch using low_rate_pitch_enc here,
        % but keeping raw index is safer for now if table missing.
    end

    % =====================================================================
    % 4. BPVC & Jitter
    % =====================================================================
    % Jitter
    if voicing_flags(1) == 0 && voicing_flags(2) == 1
        jitter_idx = 1; 
    else
        jitter_idx = 0; 
    end
    
    % BPVC
    if isempty(bpvc_input_vec)
         bpvc_idx = 0;
    else
        bpvc_val = 0;
        for i = 1:min(5, length(bpvc_input_vec))
            if bpvc_input_vec(i) > 0.5 
                bpvc_val = bitset(bpvc_val, i); 
            end
        end
        bpvc_idx = bpvc_val;
    end
end