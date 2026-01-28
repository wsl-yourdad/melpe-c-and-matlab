function best_index = gain_vq_search(target_gain_vec, gain_codebook)
% gain_vq_search: 专为6维增益矢量量化设计的搜索函数
%
% 输入:
%   target_gain_vec: (6x1) 的目标增益列矢量 (单位: dB)
%   gain_codebook:   (1024x6) 的增益码本矩阵 (单位: dB)
%
% 输出:
%   best_index: 最佳匹配码字的索引 (1 到 1024)

    % 1. 获取码本中的码字数量 (即行数)
    num_codewords = size(gain_codebook, 1);
    
    % 2. 初始化最小距离和最佳索引
    min_dist_sq = inf;
    best_index = 1;

    % 3. 遍历码本的每一行
    for i = 1:num_codewords
        
        % 从码本中取出一行，作为当前的6维码字 (1x6 的行矢量)
        current_codeword = gain_codebook(i, :);
        
        % 计算目标矢量与当前码字之间的欧氏距离的平方
        % target_gain_vec' 将 6x1 转为 1x6，以便与 codeword 进行矢量减法
        diff = target_gain_vec' - current_codeword;
        dist_sq = sum(diff.^2);
        
        % 4. 比较并更新最小距离
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq;
            best_index = i;
        end
    end

end
