function stable_lsf = stabilize_lsf(lsf, min_sep_hz, FS)
% STABILIZE_LSF 确保LSF频率严格递增且保持最小间距
%
% 输入:
%   lsf:        一个包含LSF频率的列向量 (Hz)
%   min_sep_hz: 要求的最小频率间隔 (Hz)，通常为50Hz
%   FS:         采样率

    LPC_ORD = length(lsf);
    
    % 确保LSF与边界有最小间隔
    if lsf(1) < min_sep_hz
        lsf(1) = min_sep_hz;
    end
    if lsf(LPC_ORD) > (FS/2 - min_sep_hz)
        lsf(LPC_ORD) = FS/2 - min_sep_hz;
    end
    
    % 从头到尾检查并强制拉开间距
    for i = 1:(LPC_ORD - 1)
        if lsf(i+1) < lsf(i) + min_sep_hz
            lsf(i+1) = lsf(i) + min_sep_hz;
        end
    end
    
    stable_lsf = lsf;
end
