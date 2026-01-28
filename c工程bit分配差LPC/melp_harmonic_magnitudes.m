function [fsmag_norm, fsmag_raw] = melp_harmonic_magnitudes(lpc_coeffs, pitch, fs)
% MELP_HARMONIC_MAGNITUDES 计算前10个谐波幅值并归一化
% 对应论文图 2.14 和公式 (2-17)
%
% 输入:
%   lpc_coeffs: LPC 系数向量 [1, a1, a2, ..., a10]
%   pitch:      基音周期 (采样点数) 或 基音频率 (Hz)，这里假设输入是点数(如53)
%   fs:         采样率 (通常8000)
% 输出:
%   fsmag_norm: 归一化后的10个谐波幅值 (用于量化/传输)
%   fsmag_raw:  未归一化的幅值 (用于调试)

    N_HARMONICS = 10;
    FFT_LEN = 512;
    IMPULSE_LEN = 200; % 脉冲响应长度，通常取200点(25ms)覆盖主要能量

    % --- 步骤 1: 生成合成滤波器脉冲响应 ---
    % H(z) = 1 / A(z)
    % 激励信号是一个单位脉冲 [1, 0, 0, ...]
    impulse_input = [1; zeros(IMPULSE_LEN-1, 1)];
    h = filter(1, lpc_coeffs, impulse_input);

    % --- 步骤 2: 加汉明窗 (Hamming Window) ---
    % 对应图2.14 "加窗(汉明窗)"
    w_ham = hamming(IMPULSE_LEN);
    h_windowed = h .* w_ham;

    % --- 步骤 3: 512点 FFT ---
    % 对应图2.14 "FFT(512点)"
    % 自动补零到512点
    H_fft = fft(h_windowed, FFT_LEN);
    
    % 取前一半 (0 ~ fs/2)，并取模值
    mag_spectrum = abs(H_fft(1 : FFT_LEN/2 + 1));
    
    % 频率轴 (用于查找谐波)
    freq_axis = linspace(0, fs/2, length(mag_spectrum));

    % --- 步骤 4: 峰值检测 (提取前10个谐波) ---
    % 对应图2.14 "频谱峰值检测" -> "10个最大谐波幅度值"
    
    % 计算基频 F0
    if pitch > 0
        f0 = fs / pitch; 
    else
        f0 = 0; % 清音或无效
    end
    
    fsmag_raw = ones(N_HARMONICS, 1); % 默认值为1 (平坦)

    if f0 > 60 % 只有在基频有效时才搜索 (例如 > 60Hz)
        for k = 1:N_HARMONICS
            target_freq = k * f0;
            
            if target_freq < (fs/2)
                % 1. 找到理论频率对应的 FFT 索引
                [~, center_idx] = min(abs(freq_axis - target_freq));
                
                % 2. (优化) 在理论位置附近小范围搜索局部最大值
                % 防止因为栅栏效应漏掉峰值
                search_range = 3; % 左右各搜3个点
                idx_start = max(1, center_idx - search_range);
                idx_end = min(length(mag_spectrum), center_idx + search_range);
                
                % 提取该范围内的最大值作为谐波幅值
                fsmag_raw(k) = max(mag_spectrum(idx_start:idx_end));
            else
                % 超过奈奎斯特频率，置为0或保持上一值，通常置很小的值
                fsmag_raw(k) = 0.0001; 
            end
        end
    end

    % --- 步骤 5: 归一化处理 (公式 2-17) ---
    % 公式: fsmag_i = sqrt( N / (avg + 0.0001) ) * fsmag_i
    % 这里的关键是 avg 的定义。为了让RMS归一化为1，通常 avg 指的是"平方和"或"均方值"。
    % 结合公式形式，如果 avg 是 "平方和" (Sum of Squares)，则结果是单位RMS向量。
    
    avg_val = sum(fsmag_raw .^ 2); % 计算平方和 (Total Energy of Harmonics)
    
    % 计算归一化因子
    % 解释: 如果 avg 是总能量，sqrt(N/avg) * vec 会让结果的 RMS = 1
    norm_factor = sqrt(N_HARMONICS / (avg_val + 0.0001));
    
    fsmag_norm = norm_factor * fsmag_raw;

end