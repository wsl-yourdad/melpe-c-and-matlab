function [indices, lsf_quant_out] = quantize_lsf_c_style(sf_lsfs, sf_v_dec, cb, prev_lsf)
% =========================================================================
%  MELPe 1200bps LSF Quantization (High Fidelity Version)
%  修复记录:
%    1. [Mode 1] 激活 4-bit 插值索引搜索 (inpCoef)，不再用占位符。
%    2. [Mode 1] 实现 C 标准的 20维加权预测公式，不再用简单平均。
%    3. [Mode 0] 实现 3帧独立 MSVQ (8+6 bits/frame)，填满 42bits 预算。
% =========================================================================

    % 1. 模式判定 (C标准 UV Config)
    uv_code = sprintf('%d%d%d', sf_v_dec(1), sf_v_dec(2), sf_v_dec(3));

    % [新增] 准备均值向量 (必须加上!)
    % 假设 cb.msvq_mean 已经在 build_final_codebooks.m 中加载并转为 Float (0~1)
    % 这里需要转回 Q15 以匹配 sf_lsfs 的量级
    % 如果你的 cb.msvq_mean 已经是 Q15，则直接用；如果是 Float，乘以 32768
    % 你的 sf_lsfs 明显是 0~32768 的 Q15 格式
    
    % 1. 备份原始 Q15 输入 (用于计算权重，因为 melp_lsf_weight 需要 Q15)
    sf_lsfs_q15 = sf_lsfs; 
    
    % 2. 归一化工作副本 (Work Copy)
    % 所有的减法、VQ搜索都用这个 lsf_work (Float)
    if max(abs(sf_lsfs(:))) > 100
        lsf_work = sf_lsfs / 32768.0; 
    else
        lsf_work = sf_lsfs; % 已经是 Float 了 (防御性编程)
    end
    
    % 3. 处理 Prev LSF (状态记忆)
    % 主循环里的 last_lsf_quant 可能是 Q15 (被 clamp 过)，这里要转回来
    if max(abs(prev_lsf(:))) > 100
        prev_lsf_float = prev_lsf / 32768.0;
    else
        prev_lsf_float = prev_lsf;
    end
    mean_lsf = cb.msvq_mean;
    if strcmp(uv_code, '001') 
        % =================================================================
        % Mode 1: 插值预测编码 (Interpolation Mode)
        % 结构: Base(9b) + Int_Idx(4b) + Res(8+6+6+6=26b) + Prot(3b) = 42b
        % =================================================================
        
        lsf_f2 = lsf_work(:, 2);
        lsf_f3 = lsf_work(:, 3);
        
        % A. 搜索基准 LSF (Base LSF) - 9 bits
        % 使用 Q15 的 lsp_uv_9 码本对第 3 帧进行粗量化
        [idx_base, lsf_base] = vq_search_simple(lsf_f3, cb.lsp_uv_9);
        
        % B. [关键修复] 搜索最优插值系数 (4 bits)
        % 遍历 inpCoef 的 16 组权重，找到预测误差最小的那一组
        best_err = inf;
        best_k = 1;
        best_pred_vec = zeros(20, 1);
        best_res_target = zeros(20, 1);
        
        % 准备权重用于误差评估
        w2 = melp_lsf_weight(lsf_f2);
        w3 = melp_lsf_weight(lsf_f3);
        weights_20d = [w2; w3];
        
        for k = 1:16
            % 从码本取出第 k 组权重 (20维: 前10维给F2, 后10维给F3)
            w_int = cb.inpCoef(:, k); 
            
            % [关键修复] C标准预测公式
            % Pred = w * Base + (1-w) * Prev
            % 注意: prev_lsf 和 lsf_base 都是 10维，w_int 是 20维
            % 我们需要把 Prev 和 Base 扩展对应到 F2 和 F3 的位置
            
            % F2 预测
            w_f2 = w_int(1:10);
            pred_f2 = w_f2 .* lsf_base + (1 - w_f2) .* prev_lsf;
            
            % F3 预测
            w_f3 = w_int(11:20);
            pred_f3 = w_f3 .* lsf_base + (1 - w_f3) .* prev_lsf;
            
            % 组合预测向量
            curr_pred = [pred_f2; pred_f3];
            
            % 计算残差
            curr_res = [lsf_f2; lsf_f3] - curr_pred;
            
            % 计算加权误差能量 (不需要完整VQ，只看未量化残差能量)
            err = sum(weights_20d .* (curr_res .^ 2));
            
            if err < best_err
                best_err = err;
                best_k = k;
                best_pred_vec = curr_pred;
                best_res_target = curr_res;
            end
        end
        
        % 锁定最佳插值索引 (0-based)
        idx_2nd = best_k - 1; 
        
        % C. 执行 20维残差量化 (使用选定的最佳残差)
        [idx_msvq, total_quant_res] = run_20d_msvq(best_res_target, weights_20d, cb);
        
        % D.
        % 重构与输出&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        % Mode 1 结尾:
         quant_combined = best_pred_vec + total_quant_res;
         lsf_quant_out = lsf_clamp_q15(quant_combined(11:20) * 32768.0); % 转回 Q15 钳位
        
        % 填充索引 (Mode 1 格式)
        indices.mode = 1;
        indices.idx1 = idx_base - 1;   % Base (9b)
        indices.idx2 = idx_2nd;        % Int (4b)
        indices.idx3 = idx_msvq.idx1;  % S1 (8b)
        indices.idx4 = idx_msvq.idx2;  % S2 (6b)
        indices.idx5 = idx_msvq.idx3;  % S3 (6b)
        indices.idx6 = idx_msvq.idx4;  % S4 (6b)
        
    else
        % =================================================================
        % Mode 0: 普通/全浊音模式 (Voiced/Unvoiced Mode)
        % [关键修复] 不再填 0，而是执行 3 帧独立量化
        % 结构: F1(8+6) + F2(8+6) + F3(8+6) = 42 bits
        % =================================================================
        
        indices.mode = 0;
        
        % --- Frame 1 ---
        w1 = melp_lsf_weight(lsf_work(:,1));
        % [修复] 目标向量必须减去均值！
        target_f1 = lsf_work(:,1) - mean_lsf; 
        [idx_f1_s1, r1] = vq_search_weighted(target_f1, cb.lsp_v_cb{1}, w1);
        [idx_f1_s2, r2] = vq_search_weighted(target_f1-r1, cb.lsp_v_cb{2}, w1);
        
        % --- Frame 2 ---
        w2 = melp_lsf_weight(lsf_work(:,2));
        target_f2 = lsf_work(:,2) - mean_lsf; % [修复]
        [idx_f2_s1, r1] = vq_search_weighted(target_f2, cb.lsp_v_cb{1}, w2);
        [idx_f2_s2, ~]  = vq_search_weighted(target_f2-r1, cb.lsp_v_cb{2}, w2);
        
        % --- Frame 3 ---
        w3 = melp_lsf_weight(lsf_work(:,3));
        target_f3 = lsf_work(:,3) - mean_lsf; % [修复]
        [idx_f3_s1, r1] = vq_search_weighted(target_f3, cb.lsp_v_cb{1}, w3);
        [idx_f3_s2, r2] = vq_search_weighted(target_f3-r1, cb.lsp_v_cb{2}, w3);
        
        % [修复] 输出预测值时，要把均值加回来 (虽然这里只影响预测器记忆)
        % Mode 0 结尾:
        lsf_float_out = mean_lsf + r1 + r2;
        lsf_quant_out = lsf_clamp_q15(lsf_float_out * 32768.0); % 转回 Q15 钳位%**************************************************************************
        
        % 填充索引 (Mode 0 格式: 复用 idx1~idx6 字段)
        % 为了让 pack_bits 能够通用，我们这里需要约定好含义
        % 这里直接依次填入 F1, F2, F3 的索引
        indices.idx1 = idx_f1_s1 - 1; % F1 S1 (8b)
        indices.idx2 = idx_f1_s2 - 1; % F1 S2 (6b) -> 对应 pack_bits 的 "4b" 位置? 
        % ⚠️ 注意: 这里的比特位宽 (8+6) 和 Mode 1 (9+4) 不一样！
        % 在 pack_bits 中需要根据 Mode 加以区分。
        % 暂时先存下来，下一步我们去修 pack_bits。
        indices.idx3 = idx_f2_s1 - 1; % F2 S1 (8b)
        indices.idx4 = idx_f2_s2 - 1; % F2 S2 (6b)
        indices.idx5 = idx_f3_s1 - 1; % F3 S1 (8b)
        indices.idx6 = idx_f3_s2 - 1; % F3 S2 (6b)
        
    end
end

% --- 子函数: 执行 20维 MSVQ ---
function [indices, total_res] = run_20d_msvq(target, weights, cb)
    % Stage 1 (8 bits)
    [idx1, c1] = vq_search_weighted(target, cb.res_cb{1}, weights);
    
    % Stage 2 (6 bits)
    [idx2, c2] = vq_search_weighted(target - c1, cb.res_cb{2}, weights);
    
    % Stage 3 (6 bits)
    [idx3, c3] = vq_search_weighted(target - c1 - c2, cb.res_cb{3}, weights);
    
    % Stage 4 (6 bits)
    [idx4, c4] = vq_search_weighted(target - c1 - c2 - c3, cb.res_cb{4}, weights);
    
    total_res = c1 + c2 + c3 + c4;
    indices.idx1 = idx1 - 1;
    indices.idx2 = idx2 - 1;
    indices.idx3 = idx3 - 1;
    indices.idx4 = idx4 - 1;
end

% --- 辅助函数: 简单 VQ ---
function [best_idx, best_vec] = vq_search_simple(target, codebook)
    [~, N] = size(codebook);
    min_dist = inf; best_idx = 1;
    for i=1:N
        d = sum((target - codebook(:,i)).^2);
        if d < min_dist, min_dist = d; best_idx = i; end
    end
    best_vec = codebook(:, best_idx);
end

% --- 辅助函数: 加权 VQ ---
function [best_idx, best_vec] = vq_search_weighted(target, codebook, weights)
    [~, N] = size(codebook);
    min_dist = inf; best_idx = 1;
    for i=1:N
        % 加权欧氏距离
        d = sum(weights .* (target - codebook(:,i)).^2);
        if d < min_dist, min_dist = d; best_idx = i; end
    end
    best_vec = codebook(:, best_idx);
end

% --- 辅助函数: 权重计算 ---
function w = melp_lsf_weight(lsf)
    w = ones(10,1);
    lsf = sort(lsf);
    for i=1:10
        if i==1, d1 = lsf(i); else, d1 = lsf(i)-lsf(i-1); end
        if i==10, d2 = 32767-lsf(i); else, d2 = lsf(i+1)-lsf(i); end
        if d1<50, d1=50; end
        if d2<50, d2=50; end
        w(i) = 1.0/sqrt(d1) + 1.0/sqrt(d2);
    end
end

% --- 辅助函数: 钳位 ---
function lsf = lsf_clamp_q15(lsf)
    lsf = sort(lsf);
    lsf(lsf<50) = 50;
    lsf(lsf>32700) = 32700;
    for i=1:9
        if lsf(i+1)-lsf(i) < 50
             mid = (lsf(i+1)+lsf(i))/2;
             lsf(i) = mid - 25; lsf(i+1) = mid + 25;
        end
    end
end