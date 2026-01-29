function lsf = lsf_stabilize(lsf)
    lsf = sort(lsf(:));
    min_dist = 60; % 约 7Hz 最小间距
    for k = 1:9
        if (lsf(k+1) - lsf(k)) < min_dist
            lsf(k+1) = lsf(k) + min_dist;
        end
    end
    lsf(lsf < 100) = 100; lsf(lsf > 32600) = 32600;
end