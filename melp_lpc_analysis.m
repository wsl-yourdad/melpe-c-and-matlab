function [a_i, res_frame] = melp_lpc_analysis(frame, p)
    % 执行LPC分析，获取10个预测系数 a_i 以及残差信号 res_frame
    %
    % 输入:
    %   frame: 当前的语音帧 (例如 180 点)
    %   p:     LPC的阶数 (在MELP中为 10)
    %
    % 输出:
    %   a_i:       1x10 的行向量，包含了10个LPC预测系数
    %   res_frame: 1xN 的残差信号 (由逆滤波器产生)

    % 确保 frame 是列向量
    if size(frame, 2) > 1
        frame = frame';
    end
    
    frame_len = length(frame);
    
    % 1. 应用Hamming窗 (减少边缘效应)
    windowed_frame = frame .* hamming(frame_len);

    % 2. 使用MATLAB内置函数计算LPC系数
    % a_z 的格式为 [1, a(1), a(2), ..., a(p)]
    [a_z, ~] = lpc(windowed_frame, p);

    % 3. 提取预测系数 (用于量化)
    % 根据您的 main 函数习惯，取 a(1)~a(p)
    a_i = a_z(2:end); 

    % 4. 计算残差信号 [关键补全]
    % 残差信号是通过 FIR 滤波器 A(z) = 1 - sum(a_i * z^-i) 得到的 [cite: 560]
    % 在 MATLAB 中，filter(a_z, 1, frame) 即可直接实现该逆滤波器
    res_frame = filter(a_z, 1, frame); 
end