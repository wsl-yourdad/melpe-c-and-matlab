function [indices, lsf_quant_norm] = melp_lsf_quantize_1200(lsf_hz, codebooks, fs)
    % 1. 转换为归一化频率 (0 ~ 1.0)
    % 必须除以 32768 以匹配 C 源码中的整数映射量级
    lsf_norm = lsf_hz(:) / 32768; 

    % 2. 准备权重 (10x1)
    weights = ones(10, 1);

    % 3. 执行 MSVQ 搜索 
    [indices, lsf_quant_norm] = melp_msvq_search(lsf_norm, codebooks.msvq_cb, weights);
    
    lsf_quant_norm = lsf_quant_norm(:); 
end