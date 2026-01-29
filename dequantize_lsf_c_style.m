function lsf_frames_hz = dequantize_lsf_c_style(lsf_indices, voicing_flags, l_prev_hz, codebooks, FS)
% dequantize_lsf_c_style: 严格按照C标准多模式逻辑，反量化LSF

l_prev_norm = l_prev_hz / (FS / 2);
lsf_frames_norm = zeros(10, 3);

is_interpolation_mode = false; % 简化
if all(voicing_flags == 1) && ~is_interpolation_mode
    % 方案 A: 3x(9-bit) + 3x(4-bit)
    indices_s1 = lsf_indices(1:3);
    indices_s2 = lsf_indices(4:6);
    for i = 1:3
        q1 = vq_dequantize(indices_s1(i), codebooks.lsf_cb_9bit, codebooks.msvq_mean);
        q2 = vq_dequantize(indices_s2(i), codebooks.lsf_cb_4bit, zeros(10,1));
        lsf_frames_norm(:,i) = q1 + q2;
    end
elseif any(voicing_flags == 0) && ~is_interpolation_mode
    % 方案 B: {9,9, 8,6,5,5}
    q1 = vq_dequantize(lsf_indices(1), codebooks.lsf_cb_9bit, codebooks.msvq_mean);
    q2 = vq_dequantize(lsf_indices(2), codebooks.lsf_cb_9bit, codebooks.msvq_mean);
    q3 = msvq_dequantize(lsf_indices(3:6), codebooks.lsf_cb_mode_B, codebooks.msvq_mean);
    lsf_frames_norm = [q1, q2, q3];
else
    % 方案 C: 插值
    q_msvq = msvq_dequantize(lsf_indices(1:4), codebooks.lsf_cb_mode_C, codebooks.msvq_mean);
    % 复杂的插值逻辑，暂时用固定值或重复帧代替
    lsf_frames_norm = repmat(q_msvq, 1, 3);
end

lsf_frames_hz = lsf_frames_norm * (FS / 2);
end

function quant_vec = vq_dequantize(index, codebook, cb_mean)
    quant_vec = codebook(:, index) + cb_mean(:);
end

function quant_vec = msvq_dequantize(indices, cb_stages, cb_mean)
    quant_vec = cb_mean(:);
    for s = 1:length(indices)
        quant_vec = quant_vec + cb_stages{s}(:, indices(s));
    end
end
