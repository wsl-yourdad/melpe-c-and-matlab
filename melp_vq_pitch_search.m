function best_idx = melp_vq_pitch_search(target_vec, codebook)
% 功能：针对矢量量化的最小平方误差 (MSE) 搜索 (1200bps 终极兼容版)

    % 1. 自动处理 Cell 包装 (如 res_vq{1})
    if iscell(codebook)
        codebook = codebook{1};
    end

    % 2. 强制目标为列向量 (10x1)
    target_vec = target_vec(:); 
    dim = length(target_vec);

    % 3. 【核心修正】自动对齐码本维度
    % 如果码本行数 (size1) 不等于目标维度 (10)，但列数 (size2) 等于 10，
    % 说明码本是 [256 x 10] 的行存储格式，必须转置为 [10 x 256]
    if size(codebook, 1) ~= dim && size(codebook, 2) == dim
        codebook = codebook'; 
    end

    % 4. 获取码字数量 (现在码本必然是 dim x N)
    num_codewords = size(codebook, 2);
    min_dist = inf;
    best_idx = 1;

    % 5. 遍历搜索
    for i = 1:num_codewords
        % 此时 target_vec 和 codebook(:, i) 均为 dim x 1 列向量

        % % --- 调试：仅在第一个码字时输出一次 ---
        % if i == 1
        %     fprintf('  函数内部: target_vec 维度 %s, 码本列 1 维度 %s\n', ...
        %             mat2str(size(target_vec)), mat2str(size(codebook(:, 1))));
        % end
        % % ------------------------------------

        diff = target_vec - codebook(:, i);
        dist = sum(diff.^2);
        
        if dist < min_dist
            min_dist = dist;
            best_idx = i;
        end
    end
end