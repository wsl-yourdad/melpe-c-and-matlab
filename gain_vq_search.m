function best_index = gain_vq_search(target_gain_vec, gain_codebook)
% gain_vq_search: 6维增益矢量量化 (C-Style Q8)
%
% 输入:
%   target_gain_vec: (6x1) dB Float
%   gain_codebook:   (1024x6) dB Float (loaded as Q8/256)

    % 1. 模拟 C 的 Q8 转换 (dB -> Q8 Integer)
    % 乘以 256 并四舍五入
    target_q8 = round(target_gain_vec * 256.0);
    
    % 假设 gain_codebook 载入时除以了 256，现在还原回 Q8 整数进行比较
    % 这样能保证距离计算完全一致
    cb_q8 = round(gain_codebook * 256.0); 

    % 2. 搜索
    num_codewords = size(cb_q8, 1);
    min_dist_sq = inf;
    best_index = 1;

    for i = 1:num_codewords
        % 取出整数码字
        current_codeword = cb_q8(i, :);
        
        % 整数域减法
        diff = target_q8' - current_codeword;
        
        % 整数平方和 (模拟 C 的 mult/mac)
        dist_sq = sum(diff.^2);
        
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq;
            best_index = i;
        end
    end

end