function [indices, lsf_quant_out] = quantize_lsf_c_style(sf_lsfs, sf_v_dec, cb, prev_lsf)
% =========================================================================
%  MELPe 1200bps LSF Quantization (C-Style Bit-Exact Version)
%  修正:
%    1. 修正 Mode 1 触发条件: 包含 000(VVV), 001, 010, 100 等混合模式
%    2. 引入 C 风格 Q11 权重计算 (simulated in float)
%    3. 严格对齐插值与残差搜索逻辑
% =========================================================================

    % 1. 模式判定 (匹配 C 代码 switch-case 逻辑)
    % C代码 uv_config 位序: MSB=F1, LSB=F3 (即 uv_flag[0]...uv_flag[2])
    % 假设 sf_v_dec = [uv1, uv2, uv3] (1=Unvoiced, 0=Voiced in C? Wait.)
    % ⚠️ 注意: C 代码中 par[i].uv_flag=1 表示清音(Unvoiced)。
    %    MATLAB 中 sf_v_dec=1 表示浊音(Voiced), 0 表示清音(Unvoiced)?
    %    需确认您的 sf_v_dec 定义。通常 MELP MATLAB 仿真里 1=Voiced。
    %    如果 sf_v_dec: 1=Voiced, 0=Unvoiced。
    %    则 C 的 uv_flag = 1 - sf_v_dec。
    
    % 转换逻辑: 构造 C 风格 uv_config (3 bits)
    % C: 1=UV, 0=V. 
    c_uv_flags = 1 - sf_v_dec; % 转换: 1->0(V), 0->1(UV)
    uv_config = c_uv_flags(1)*4 + c_uv_flags(2)*2 + c_uv_flags(3)*1;

    % 2. 准备工作数据
    % 假设输入 sf_lsfs 是 Q15 (0~32768)
    sf_lsfs_q15 = sf_lsfs;
    
    % 转为 Float (0~1) 用于计算，匹配码本量级
    lsf_float = sf_lsfs / 32768.0; 
    prev_lsf_float = prev_lsf;
    if max(prev_lsf_float) > 100, prev_lsf_float = prev_lsf_float / 32768.0; end
    
    mean_lsf = cb.msvq_mean; % Float 0~1

    % =====================================================================
    % 逻辑分流 (参考 qnt12.c lsf_vq)
    % =====================================================================
    % Case 7 (111): All UV
    % Case 6 (110): UV UV V
    % Case 5 (101): UV V UV
    % Case 3 (011): V UV UV
    % Default (000, 001, 010, 100): Mode 1 (Interpolation)
    
    use_mode_1 = false;
    if ismember(uv_config, [0, 1, 2, 4]) 
        use_mode_1 = true; % VVV, VVU, VUV, UVV -> 插值模式
    end

if use_mode_1
        % =================================================================
        % Mode 1: 插值预测编码 (42 bits) - C-Style Bit-Exact 逻辑
        % =================================================================
        indices.mode = 1;
        lsf_f2 = lsf_float(:, 2);
        lsf_f3 = lsf_float(:, 3);
        
        % --- 分支逻辑: 001 (Mixed) vs Others (VVV) ---
        % 注意: uv_config 定义 (C代码): 1=001(VVU). 其他通常是 VVV.
        
        if uv_config == 1  % Case: 001 (VVU)
            is_case_001 = true;
            
            % A. Base LSF (9 bits)
            [idx_base, lsf_base] = vq_search_simple(lsf_f3, cb.lsp_uv_9);
            
            % 记录 Base 索引
            indices.idx1 = idx_base - 1; 
            
        else % Case: VVV (000), 010, 100
            is_case_001 = false;
            
            % A. Base LSF (24 bits = 8+6+5+5)
            % 使用 lsp_v_cb (2400bps 码本)
            target_base = lsf_f3 - mean_lsf;
            w3 = melp_lsf_weight_c_style(sf_lsfs_q15(:,3));
            
            [idx_b1, r1] = vq_search_weighted(target_base, cb.lsp_v_cb{1}, w3);
            [idx_b2, r2] = vq_search_weighted(target_base-r1, cb.lsp_v_cb{2}, w3);
            [idx_b3, r3] = vq_search_weighted(target_base-r1-r2, cb.lsp_v_cb{3}, w3);
            [idx_b4, r4] = vq_search_weighted(target_base-r1-r2-r3, cb.lsp_v_cb{4}, w3);
            
            lsf_base = mean_lsf + r1 + r2 + r3 + r4;
            
            % 记录 Base 索引 (复用 idx1 字段结构，打包时拆分)
            indices.base_stages = [idx_b1-1, idx_b2-1, idx_b3-1, idx_b4-1];
        end
        
        % B. Search Interpolation Weights (4 bits)
        % -------------------------------------------------------------
        % (这部分逻辑通用，但 Residual 量化级数不同)
        w2 = melp_lsf_weight_c_style(sf_lsfs_q15(:,2));
        w3 = melp_lsf_weight_c_style(sf_lsfs_q15(:,3));
        weights_20d = [w2; w3];

        best_err = inf;
        best_k = 1;
        best_pred_vec = zeros(20, 1);
        best_res_target = zeros(20, 1);
        
        for k = 1:16
            w_int = cb.inpCoef(:, k);
            w_f2 = w_int(1:10);
            w_f3 = w_int(11:20);
            
            pred_f2 = w_f2 .* lsf_base + (1 - w_f2) .* prev_lsf_float;
            pred_f3 = w_f3 .* lsf_base + (1 - w_f3) .* prev_lsf_float;
            
            curr_pred = [pred_f2; pred_f3];
            curr_res = [lsf_f2; lsf_f3] - curr_pred; 
            
            err = sum(weights_20d .* (curr_res .^ 2));
            if err < best_err
                best_err = err;
                best_k = k;
                best_pred_vec = curr_pred;
                best_res_target = curr_res;
            end
        end
        indices.idx2 = best_k - 1; % Interpolation Index

        % C. Residual Quantization (差异化)
        % -------------------------------------------------------------
        if is_case_001
            % 001 模式: 4级残差 (26 bits)
            [idx_msvq, total_quant_res] = run_20d_msvq(best_res_target, weights_20d, cb, 4);
        else
            % VVV 模式: 仅2级残差 (14 bits)
            [idx_msvq, total_quant_res] = run_20d_msvq(best_res_target, weights_20d, cb, 2);
        end
        
        indices.idx3 = idx_msvq.idx1;
        indices.idx4 = idx_msvq.idx2;
        indices.idx5 = idx_msvq.idx3; % VVV模式下为0
        indices.idx6 = idx_msvq.idx4; % VVV模式下为0

        % Output Construction
        quant_combined = best_pred_vec + total_quant_res;
        lsf_quant_out = lsf_clamp_q15(quant_combined(11:20) * 32768.0);
        
        % 标记模式以便 Packer 识别
        indices.is_001 = is_case_001;

    else
        % =================================================================
        % Mode 0: 独立量化 (42 bits)
        % 适用于: UV UV UV (111), UV UV V (110) 等
        % =================================================================
        indices.mode = 0;
        
        % Frame 1
        w1 = melp_lsf_weight_c_style(sf_lsfs_q15(:,1));
        target_f1 = lsf_float(:,1) - mean_lsf;
        [idx_f1_s1, r1] = vq_search_weighted(target_f1, cb.lsp_v_cb{1}, w1);
        [idx_f1_s2, r2] = vq_search_weighted(target_f1-r1, cb.lsp_v_cb{2}, w1);
        
        % Frame 2
        w2 = melp_lsf_weight_c_style(sf_lsfs_q15(:,2));
        target_f2 = lsf_float(:,2) - mean_lsf;
        [idx_f2_s1, r1] = vq_search_weighted(target_f2, cb.lsp_v_cb{1}, w2);
        [idx_f2_s2, ~]  = vq_search_weighted(target_f2-r1, cb.lsp_v_cb{2}, w2);
        
        % Frame 3
        w3 = melp_lsf_weight_c_style(sf_lsfs_q15(:,3));
        target_f3 = lsf_float(:,3) - mean_lsf;
        [idx_f3_s1, r1] = vq_search_weighted(target_f3, cb.lsp_v_cb{1}, w3);
        [idx_f3_s2, r2] = vq_search_weighted(target_f3-r1, cb.lsp_v_cb{2}, w3);
        
        lsf_float_out = mean_lsf + r1 + r2;
        lsf_quant_out = lsf_clamp_q15(lsf_float_out * 32768.0);
        
        indices.idx1 = idx_f1_s1 - 1;
        indices.idx2 = idx_f1_s2 - 1;
        indices.idx3 = idx_f2_s1 - 1;
        indices.idx4 = idx_f2_s2 - 1;
        indices.idx5 = idx_f3_s1 - 1;
        indices.idx6 = idx_f3_s2 - 1;
    end
end

% --- C-Style LSF Weighting (Q11 scale simulated in Float) ---
function w = melp_lsf_weight_c_style(lsf_q15)
    % Input: lsf_q15 (10x1 Integer 0~32768)
    % Output: w (10x1 Float, scaled to match C's Q11 magnitude approx)
    
    lsf = sort(lsf_q15);
    w = zeros(10,1);
    
    % C Code: 1/sqrt(d)
    % C Q11: 2048.0
    
    for i=1:10
        if i==1, d1 = lsf(i); else, d1 = lsf(i)-lsf(i-1); end
        if i==10, d2 = 32767-lsf(i); else, d2 = lsf(i+1)-lsf(i); end
        
        if d1<50, d1=50; end
        if d2<50, d2=50; end
        
        % MATLAB Float simulation of C logic
        % w[i] = (1/sqrt(d1) + 1/sqrt(d2)) * Scaling
        % In C, weights are relative. 
        % Let's use standard formula but ensure large dynamic range isn't clamped
        val = 1.0/sqrt(double(d1)) + 1.0/sqrt(double(d2));
        
        % Scale up to avoid underflow in distance calc, mimicking Q11 (x2048)
        w(i) = val * 2048.0; 
    end
end

% --- MSVQ Search (Unchanged but ensuring weights used) ---
function [indices, total_res] = run_20d_msvq(target, weights, cb, stages)
    % 支持指定级数 (stages = 2 或 4)
    [idx1, c1] = vq_search_weighted(target, cb.res_cb{1}, weights);
    [idx2, c2] = vq_search_weighted(target - c1, cb.res_cb{2}, weights);
    
    if stages >= 3
        [idx3, c3] = vq_search_weighted(target - c1 - c2, cb.res_cb{3}, weights);
    else
        idx3 = 1; c3 = 0;
    end
    
    if stages >= 4
        [idx4, c4] = vq_search_weighted(target - c1 - c2 - c3, cb.res_cb{4}, weights);
    else
        idx4 = 1; c4 = 0;
    end
    
    total_res = c1 + c2 + c3 + c4;
    indices.idx1 = idx1 - 1;
    indices.idx2 = idx2 - 1;
    indices.idx3 = idx3 - 1;
    indices.idx4 = idx4 - 1;
end

% --- Weighted VQ ---
function [best_idx, best_vec] = vq_search_weighted(target, codebook, weights)
    [~, N] = size(codebook);
    min_dist = inf; best_idx = 1;
    for i=1:N
        diff = target - codebook(:,i);
        % Weighted MSE: sum(w * diff^2)
        d = sum(weights .* (diff.^2));
        if d < min_dist, min_dist = d; best_idx = i; end
    end
    best_vec = codebook(:, best_idx);
end

function [best_idx, best_vec] = vq_search_simple(target, codebook)
    [~, N] = size(codebook);
    min_dist = inf; best_idx = 1;
    for i=1:N
        d = sum((target - codebook(:,i)).^2);
        if d < min_dist, min_dist = d; best_idx = i; end
    end
    best_vec = codebook(:, best_idx);
end

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