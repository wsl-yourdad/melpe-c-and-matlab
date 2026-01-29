function [G1_idx, G2_idx] = melp_quantize_gain(G1_unquantized, G2_unquantized, G2_previous)
    % 实现 G1 和 G2 的增益量化
    % 基于论文 和 流程图 2.13
    
    % G1_unquantized: 当前帧未量化的G1 (来自melp_gain_calculator)
    % G2_unquantized: 当前帧未量化的G2 (来自melp_gain_calculator)
    % G2_previous:    前一帧(G2p)未量化的G2

    % -------------------------------------------------
    % 步骤 1: G2 (第二增益) 量化
    % -------------------------------------------------
    % "G2用 5bit 去进行均匀量化，量化范围为 10dB 到 77dB"
    gain_min_g2 = 10;
    gain_max_g2 = 77;
    num_levels_g2 = 32; % 5 bits = 2^5 = 32 个电平 (索引 0-31)
    
    % 1.1 将 G2 钳位(Clamping)到有效范围内
    G2_clamped = max(gain_min_g2, G2_unquantized); %>=10
    G2_clamped = min(gain_max_g2, G2_clamped);     %<=77
    
    % 1.2 执行均匀量化：将 [10, 77] 范围的值 映射到 [0, 31] 的索引
    G2_idx = round( ((G2_clamped - gain_min_g2) / (gain_max_g2 - gain_min_g2)) * (num_levels_g2 - 1) );

    % -------------------------------------------------
    % 步骤 2: G1 (第一增益) 量化 (严格按图2.13)
    % -------------------------------------------------
    G2p = G2_previous;
    G2 = G2_unquantized;
    G1 = G1_unquantized;

    % 2.1 检查特殊 "索引0" 条件
    % "|G2-G2p|<5 & |G1-(G2+G2p)/2|<3"
    if (abs(G2 - G2p) < 5) && (abs(G1 - (G2 + G2p) / 2) < 3)
        % 2.1.1 Y 分支: "量化索引置0"
        G1_idx = 0;
        % (根据论文，这个索引0会告诉解码器将 G1 设为 G2 和 G2p 的均值)
        
    else
        % 2.1.2 N 分支: 执行3-bit自适应量化
        
        % 2.2 计算自适应量化范围
        % "gain_max = max(G2p, G2) + 6"
        % "gain_min = min(G2p, G2) - 6"
        gain_max = max(G2p, G2) + 6;
        gain_min = min(G2p, G2) - 6;
        
        % 2.3 钳位自适应范围
        % "gain_min < 10" -> "gain_min = 10"
        if gain_min < 10
            gain_min = 10;
        end
        
        % "gain_max > 77" -> "gain_max = 77"
        if gain_max > 77
            gain_max = 77;
        end
        
        % 2.4 "对增益进行7电平的均匀量化"
        % (这7个电平 + "索引0" = 8个总电平 = 3 bit)
        num_levels_g1_adaptive = 7; % 7 个电平 (用于索引 1-7)
        
        % 2.5 将 G1 钳位到自适应范围内
        G1_clamped = max(gain_min, G1);
        G1_clamped = min(gain_max, G1_clamped);
        
        % 2.6 执行7电平的均匀量化
        % 将 [gain_min, gain_max] 范围的值 映射到 [0, 6] 的索引
        adaptive_idx = 0; % 默认值
        if (gain_max - gain_min) > 0
            adaptive_idx = round( ((G1_clamped - gain_min) / (gain_max - gain_min)) * (num_levels_g1_adaptive - 1) );
        end
        
        % 2.7 将索引 [0, 6] 偏移到 [1, 7]
        G1_idx = adaptive_idx + 1;
    end
end