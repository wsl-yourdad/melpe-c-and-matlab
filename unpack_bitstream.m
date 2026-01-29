function p = unpack_bitstream(s)
    % 严格对应 pack_bits_c_style(pitch, gain, uv, lsf, jit, bpvc, fs)
    p.uv_index  = bin2dec(s(2:4));      % 3 bits
    p.uv_parity = bin2dec(s(5));        % 1 bit (Pad位在C代码中有时用于校验)
    p.pitch_idx = bin2dec(s(6:14));     % 9 bits
    p.lsf_bits  = s(15:56);             % 42 bits
    p.gain_idx  = bin2dec(s(57:66));    % 10 bits
    p.bpvc_bits = s(67:72);             % 6 bits (2位/帧 x 3帧)
    p.fs_idx    = bin2dec(s(73:80));    % 8 bits
    p.jitter    = bin2dec(s(81));       % 1 bit
end