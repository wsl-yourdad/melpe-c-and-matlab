function [indices, quantized_vec] = melp_msvq_search(target_vec, msvq_cell, weights)
% 功能：多级矢量量化 (MSVQ) 搜索
% 输入：
%   target_vec - [10x1] 原始 LSF 矢量
%   msvq_cell  - {1x4} 码本细胞数组
%   weights    - [10x1] 知觉加权系数 (可先设为全 1)
% 输出：
%   indices    - [1x4] 各级码本索引
%   quantized_vec - 量化还原后的矢量

    curr_residual = target_vec;
    num_stages = length(msvq_cell);
    indices = zeros(1, num_stages);
    quantized_vec = zeros(size(target_vec));

    for s = 1:num_stages

        % if s == 1
        %     fprintf('\n--- MSVQ Stage 1 Debug ---\n');
        %     fprintf('Target(1:3): [%.2f, %.2f, %.2f]\n', curr_residual(1), curr_residual(2), curr_residual(3));
        %     fprintf('CB Stage 1 Min/Max: [%.2f, %.2f]\n', min(msvq_cell{s}(:)), max(msvq_cell{s}(:)));
        % end

        cb = msvq_cell{s}; % 当前级码本 [10 x N]
        num_codewords = size(cb, 2);
        
        min_dist = inf;
        best_idx = 1;
        
        % 每一级进行加权欧氏距离搜索
        for i = 1:num_codewords
            diff = curr_residual - cb(:, i);
            dist = sum(weights .* (diff.^2));
            if dist < min_dist
                min_dist = dist;
                best_idx = i;
            end
        end
        
        indices(s) = best_idx;
        stage_vec = cb(:, best_idx);
        quantized_vec = quantized_vec + stage_vec;
        % 更新残差进入下一级
        curr_residual = curr_residual - stage_vec;
    end
end