function bit_str = pack_bits_c_style(pitch_idx, gain_idx, uv_encoded_val, lsf_indices, jit_idx, bpvc_idx, fs_idx)
% =========================================================================
%  MELPe 1200bps Bit Packer (81 bits) - 修正版
%  输入差异: 
%     uv_encoded_val: 必须是整数标量 (0-7)，不再是 [0 1 0] 向量!
%     pitch_idx: 必须是 9 bits (0-511)
% =========================================================================

    % 1. Sync Bit (1 bit)
    b_sync = '0'; 

    % 2. UV (3 bits) & Parity (1 bit)
    % 输入已经是编码好的 index (0-7)，直接转二进制
    b_uv = dec2bin(uv_encoded_val, 3);
    
    % 计算 Parity (C代码逻辑: parity of uv_index)
    % 简单奇偶校验: 1的个数是奇数则为1? 
    % MELPe C代码: parity(x,3) returns 1 if odd number of 1s.
    parity_val = mod(sum(b_uv == '1'), 2);
    b_uv_parity = dec2bin(parity_val, 1);

    % 3. Pitch (9 bits)
    % 此时 pitch_idx 已经被 quantize 函数截断为 9 bits，不会溢出
    b_pitch = dec2bin(pitch_idx, 9);

    % 4. LSF (42 bits) - 保持不变
    % =====================================================================
b_lsf_full = '';
    
    if isfield(lsf_indices, 'mode') && lsf_indices.mode == 1
        % Mode 1
        
        % 判断是 001 还是 VVV (通过之前记录的 flag)
        is_001 = isfield(lsf_indices, 'is_001') && lsf_indices.is_001;
        
        if is_001
            % Case 001: Base(9) + Int(4) + Res(8+6+6+6) + Prot(3)
            b_base = dec2bin(lsf_indices.idx1, 9);
            b_int  = dec2bin(lsf_indices.idx2, 4);
            
            b_s1   = dec2bin(lsf_indices.idx3, 8);
            b_s2   = dec2bin(lsf_indices.idx4, 6);
            b_s3   = dec2bin(lsf_indices.idx5, 6);
            b_s4   = dec2bin(lsf_indices.idx6, 6);
            b_prot = '000'; % 3 bit protection
            
            b_lsf_full = [b_base, b_int, b_s1, b_s2, b_s3, b_s4, b_prot];
        else
            % Case VVV: Base(8+6+5+5) + Int(4) + Res(8+6)
            stages = lsf_indices.base_stages; % [s1, s2, s3, s4]
            b_base_s1 = dec2bin(stages(1), 8);
            b_base_s2 = dec2bin(stages(2), 6);
            b_base_s3 = dec2bin(stages(3), 5);
            b_base_s4 = dec2bin(stages(4), 5);
            
            b_int     = dec2bin(lsf_indices.idx2, 4);
            
            b_res_s1  = dec2bin(lsf_indices.idx3, 8);
            b_res_s2  = dec2bin(lsf_indices.idx4, 6);
            
            b_lsf_full = [b_base_s1, b_base_s2, b_base_s3, b_base_s4, ...
                          b_int, ...
                          b_res_s1, b_res_s2];
        end
    else
        % Mode 0: (8+6) * 3
        b_f1_s1 = dec2bin(lsf_indices.idx1, 8);
        b_f1_s2 = dec2bin(lsf_indices.idx2, 6);
        b_f2_s1 = dec2bin(lsf_indices.idx3, 8);
        b_f2_s2 = dec2bin(lsf_indices.idx4, 6);
        b_f3_s1 = dec2bin(lsf_indices.idx5, 8);
        b_f3_s2 = dec2bin(lsf_indices.idx6, 6);
        b_lsf_full = [b_f1_s1, b_f1_s2, b_f2_s1, b_f2_s2, b_f3_s1, b_f3_s2];
    end
    
    if length(b_lsf_full) ~= 42
        error('LSF打包长度错误! 当前模式: %d, 长度: %d', lsf_indices.mode, length(b_lsf_full));
    end

    % 5. Gain (10 bits)
    b_gain = dec2bin(gain_idx, 10);

    % 6. BPVC (6 bits)
    b_bpvc = dec2bin(bpvc_idx, 6);

    % 7. FS Magnitude (8 bits)
    b_fs = dec2bin(fs_idx, 8);

    % 8. Jitter (1 bit)
    b_jit = dec2bin(jit_idx, 1);

    % --- 拼接 (81 bits) ---
    bit_str = [b_sync, b_uv, b_uv_parity, b_pitch, b_lsf_full, b_gain, b_bpvc, b_fs, b_jit];
    
    if length(bit_str) ~= 81
        warning('Packer长度异常: %d bits (预期 81)', length(bit_str));
    end
end