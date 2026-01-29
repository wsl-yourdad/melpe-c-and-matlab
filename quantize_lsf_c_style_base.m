function [lsf_indices, lsf_quant_final_out] = quantize_lsf_c_style(lsf_vectors, voicing_flags, codebooks, lsf_last)
% [lsf_indices, lsf_quant_final_out] = quantize_lsf_c_style(...)
%
% 输入:
%   lsf_vectors:   [10x3] 矩阵，必须是 Q15 整数格式 (0~32767)！
%   voicing_flags: [1x3]  向量 (0=清, 1=浊)
%   codebooks:     THE_FINAL_CODEBOOK.mat 里的结构体
%   lsf_last:      [10x1] 向量，上一帧的 Q15 量化值
%
% 输出:
%   lsf_indices:         [1xN] 索引 (1-based, 编码时需减1)
%   lsf_quant_final_out: [10x1] Q15 格式的量化值 (用于更新状态)

    % 1. 准备数据
    % 这里的输入已经是 Q15 (由主循环保证)，所以直接使用
    lsf_vec_q15 = lsf_vectors; 
    lsf_last_q15 = lsf_last;

    % 2. 获取均值 (Q15 整数)
    if isfield(codebooks, 'msvq_mean')
        cb_mean = codebooks.msvq_mean;
    else
        error('❌ 码本中缺少 msvq_mean 字段！');
    end

    % =========================================================================
    % 逻辑分支：根据 SC1200 标准决定编码路径
    % =========================================================================
    
    % 只有 "全清音" (UV-UV-UV) 走路径1，其他所有情况(包含混合)走路径2
    is_all_uv = all(voicing_flags == 0);
    
    if is_all_uv
        % === [路径 1] 全清音模式 (42 bits) ===
        % 结构: Frame1(9b) + Frame2(9b) + Frame3(24b)
        
        % 1. 量化第1帧 (9 bits) - 绝对量化
        % lsf_cb_9bit 是 [10 x 512] 的 Q15 整数矩阵
        [idx1, ~] = vq_search(lsf_vec_q15(:,1), codebooks.lsf_cb_9bit, cb_mean);
        
        % 2. 量化第2帧 (9 bits) - 绝对量化
        [idx2, ~] = vq_search(lsf_vec_q15(:,2), codebooks.lsf_cb_9bit, cb_mean);
        
        % 3. 量化第3帧 (24 bits = 8+6+5+5) - MSVQ
        % 这里的 target 是 (LSF - Mean)
        [indices_f3, lsf_q3] = msvq_search_4stage(lsf_vec_q15(:,3), codebooks.lsf_cb_mode_B, cb_mean);
        
        % 4. 打包索引
        lsf_indices = [idx1, idx2, indices_f3];
        
        % 5. 【关键修正】直接返回 Q15 整数，不要乘 FS/2！
        lsf_quant_final_out = lsf_q3; 
        
    else
        % === [路径 2] 浊音/混合模式 (42 bits) ===
        % 结构: Frame3(24b) + Interp(4b) + Residual(14b)
        
        % 1. 先量化第3帧 (End Frame) - 24 bits
        [indices_f3, lsf_q3] = msvq_search_4stage(lsf_vec_q15(:,3), codebooks.lsf_cb_mode_B, cb_mean);
        
        % 2. 插值搜索 (Interpolation Search) - 4 bits
        % 利用 lsf_last (Q15) 和 lsf_q3 (Q15) 拟合中间帧
        best_int_idx = 1;
        min_err = inf;
        best_pred_f1 = zeros(10,1);
        best_pred_f2 = zeros(10,1);
        
        % inpCoef 是 normalized float (0~1)，这是对的
        inp_matrix = codebooks.inpCoef; 
        
        for k = 1:16
            w_f1 = inp_matrix(k, 1:10)'; 
            w_f2 = inp_matrix(k, 11:20)';
            
            % 预测结果也是 Q15
            p_f1 = (1 - w_f1) .* lsf_last_q15 + w_f1 .* lsf_q3;
            p_f2 = (1 - w_f2) .* lsf_last_q15 + w_f2 .* lsf_q3;
            
            % 计算误差 (Q15域)
            err = sum((lsf_vec_q15(:,1) - p_f1).^2) + sum((lsf_vec_q15(:,2) - p_f2).^2);
            
            if err < min_err
                min_err = err;
                best_int_idx = k;
                best_pred_f1 = p_f1;
                best_pred_f2 = p_f2;
            end
        end
        
        % 3. 残差量化 (Residual Quantization) - 14 bits
        % 残差 = 目标(Q15) - 预测(Q15) = 残差(Q15)
        res_target_f1 = lsf_vec_q15(:,1) - best_pred_f1;
        res_target_f2 = lsf_vec_q15(:,2) - best_pred_f2;
        
        % Stage 1 (8 bits): Frame 1 Residual
        % 注意：残差没有均值，所以 mean 传 0
        [res_idx1, ~] = vq_search(res_target_f1, codebooks.lsf_cb_mode_B{1}, zeros(10,1));
        
        % Stage 2 (6 bits): Frame 2 Residual
        [res_idx2, ~] = vq_search(res_target_f2, codebooks.lsf_cb_mode_B{2}, zeros(10,1));
        
        % 4. 打包索引
        lsf_indices = [indices_f3, best_int_idx, res_idx1, res_idx2];
        
        % 5. 【关键修正】直接返回 Q15 整数
        lsf_quant_final_out = lsf_q3;
    end
end

% --- 内部辅助函数 1: 单级 VQ 搜索 ---
function [best_idx, quant_vec] = vq_search(target_vec, codebook, cb_mean)
    % 1. 将输入去均值 (如果是绝对量化)
    target_residual = target_vec(:) - cb_mean(:);
    
    % 2. 搜索 (Euclidean Distance)
    % codebook 列向量是 Q15，target_residual 也是 Q15
    [~, best_idx] = min(sum((codebook - target_residual).^2, 1));
    
    % 3. 重建量化值 (加上均值)
    quant_vec = codebook(:, best_idx) + cb_mean(:);
end

% --- 内部辅助函数 2: 4级 MSVQ 搜索 ---
function [indices, quantized_vec_final] = msvq_search_4stage(target_vec, cb_stages, cb_mean)
    % 1. 初始残差 = 目标 - 均值
    residual = target_vec(:) - cb_mean(:);
    indices = zeros(1, 4);
    
    quantized_vec_accum = zeros(size(residual)); % 累积量化残差
    
    num_stages = min(4, length(cb_stages));
    
    for s = 1:num_stages
        cb = cb_stages{s};
        
        % 每一级都在逼近当前的 residual
        % 这里的 cb 本身就是残差码本，不需要再减均值
        [idx, ~] = vq_search(residual, cb, zeros(10,1));
        
        indices(s) = idx;
        
        % 提取当前级的最佳向量
        stage_vec_actual = cb(:, idx);
        
        % 累加到总残差估计中
        quantized_vec_accum = quantized_vec_accum + stage_vec_actual;
        
        % 更新剩余残差
        residual = residual - stage_vec_actual;
    end
    
    % 最终重建值 = 均值 + 所有级的残差和
    quantized_vec_final = cb_mean(:) + quantized_vec_accum;
end