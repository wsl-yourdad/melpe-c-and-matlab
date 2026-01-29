function G_db = melp_calculate_gain(window)
    % 严格实现论文公式 (2-12)增益公式
    % G = 10*log10(0.01 + (1/L) * sum(S_n^2))

    L = length(window);

    % 安全检查：如果窗口为空（例如在信号开头），返回0
    if L == 0
        G_db = 0;
        return;
    end

    % (1/L) * sum(S_n^2) 是均方能量 (Mean Square Energy)
    mean_sq_energy = sum(window .^ 2) / L;

    % 应用公式 (2-12)，0.01是为了防止log(0)
    G_db = 10 * log10(0.01 + mean_sq_energy);

    % "如果增益计算出来的植小于0，则将其值设置为0"
    if G_db < 0
        G_db = 0;
    end
end

