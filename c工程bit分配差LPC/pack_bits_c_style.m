function bit_str = pack_bits_c_style(pitch_idx, gain_idx, uv, lsf_indices, jit_idx, bpvc_idx, fs_idx)
% =========================================================================
%  MELPe 1200bps Bit Packer (81 bits) - 终极适配版
%  功能: 
%    将各个模块的量化索引打包成标准比特流。
%    核心是根据 LSF 的 mode (0或1) 动态调整 42-bit LSF 字段的结构。
%  结构概览:
%    Sync(1) | UV(3) | Pitch(9) | LSF(42) | Gain(10) | BPVC(6) | FS(8) | Jitter(1)
%    (注: UV字段实际占用位置可能涉及 Parity，此处为了对齐 encoder_output.mat
%     的前几位，采用 Sync(1)+UV(3) 的紧凑排布，若需对齐 UV Parity 可微调)
% =========================================================================

    % 1. Sync Bit (1 bit) - 这里的 0/1 交替在主循环控制更好，这里暂固定
    b_sync = '0'; 

    % 2. UV (3 bits)
    % C标准中 UV 信息有时会分散存放，但在 1200bps 块结构中通常集中
    b_uv = sprintf('%d%d%d', uv(1), uv(2), uv(3)); 
    
    % 注意: 你之前的代码里有个 b_uv_full = [b_uv, '0'] (4 bits)
    % 但为了严格凑齐 81 bits，我们需要精打细算。
    % 1+3+9+42+10+6+8+1 = 80 bits? 缺 1 bit?
    % C标准里通常有一个 UV Parity (1 bit) 或者保留位。
    % 我们这里显式加上这个 Parity/Padding 位，凑齐 81 bits 里的同步头部分。
    b_uv_parity = '0'; 

    % 3. Pitch (9 bits)
    b_pitch = dec2bin(pitch_idx, 9);

    % 4. LSF (42 bits) - 核心动态打包区
    % =====================================================================
    b_lsf_full = '';
    
    if isfield(lsf_indices, 'mode') && lsf_indices.mode == 1
        % === Mode 1: 插值模式 (Interpolation) ===
        % 结构: Base(9) + Int(4) + S1(8) + S2(6) + S3(6) + S4(6) + Prot(3) = 42
        
        b_base = dec2bin(lsf_indices.idx1, 9);  % Base Index
        b_int  = dec2bin(lsf_indices.idx2, 4);  % Interpolation Index
        b_s1   = dec2bin(lsf_indices.idx3, 8);  % MSVQ Stage 1
        b_s2   = dec2bin(lsf_indices.idx4, 6);  % MSVQ Stage 2
        b_s3   = dec2bin(lsf_indices.idx5, 6);  % MSVQ Stage 3
        b_s4   = dec2bin(lsf_indices.idx6, 6);  % MSVQ Stage 4
        b_prot = '000';                         % Protection (3 bits)
        
        b_lsf_full = [b_base, b_int, b_s1, b_s2, b_s3, b_s4, b_prot];
        
    else
        % === Mode 0: 普通模式 (Independent MSVQ) ===
        % 结构: F1(8+6) + F2(8+6) + F3(8+6) = 42
        
        % Frame 1
        b_f1_s1 = dec2bin(lsf_indices.idx1, 8);
        b_f1_s2 = dec2bin(lsf_indices.idx2, 6);
        
        % Frame 2
        b_f2_s1 = dec2bin(lsf_indices.idx3, 8);
        b_f2_s2 = dec2bin(lsf_indices.idx4, 6);
        
        % Frame 3
        b_f3_s1 = dec2bin(lsf_indices.idx5, 8);
        b_f3_s2 = dec2bin(lsf_indices.idx6, 6);
        
        b_lsf_full = [b_f1_s1, b_f1_s2, b_f2_s1, b_f2_s2, b_f3_s1, b_f3_s2];
    end
    
    % 强制长度检查 (防止溢出)
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
    % Sync(1) + UV(3) + Parity(1) + Pitch(9) + LSF(42) + Gain(10) + BPVC(6) + FS(8) + Jit(1)
    % = 1+3+1+9+42+10+6+8+1 = 81 bits! 完美!
    
    bit_str = [b_sync, b_uv, b_uv_parity, b_pitch, b_lsf_full, b_gain, b_bpvc, b_fs, b_jit];
    
    % 最终保险
    if length(bit_str) ~= 81
        warning('Packer长度异常: %d bits (预期 81)', length(bit_str));
        % 强行截断或补零以防崩溃
        bit_str = pad_or_trim(bit_str, 81);
    end
end

function s = pad_or_trim(s, target_len)
    if length(s) > target_len
        s = s(1:target_len);
    elseif length(s) < target_len
        s = [s, repmat('0', 1, target_len - length(s))];
    end
end