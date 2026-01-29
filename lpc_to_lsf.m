function lsf = lpc_to_lsf(a_i, fallback_lsf)
% =========================================================================
% 功能: 将 10 阶 LPC 系数转换为缩放后的 LSF 频率值 (Scaled Hz)
% 适配: 1.2kb/s MF-MELP 仿真与 16 位定点数码本 (缩放因子 8.192)
% 
% 输入:
%   a_i          - 1x10 的 LPC 预测系数 (不含第一个系数 1) 
%   fallback_lsf - 10x1 的回退 LSF 矢量 (用于稳定性校验失败时)
%
% 输出:
%   lsf          - 10x1 的缩放 LSF 矢量 (范围 0 ~ 32768)
% =========================================================================

    persistent frame_counter;      %静态局部变量，运行结束后保留当前值
    if isempty(frame_counter), frame_counter = 0; end
    frame_counter = frame_counter + 1;

    FS = 8000;
    SCALING_FACTOR = 8.192; % 核心：将 Hz 映射到定点码本量级 (32768->4000)
    
    % 1. 默认值处理：若未传回退值，使用码本量级的典型分布值
    if nargin < 2 || isempty(fallback_lsf)
        fallback_lsf = ((1:10)' * 350) * SCALING_FACTOR; 
    end

    % 2. 构造多项式系数 [1, -a1, -a2, ..., -a10]
    % 注意：根据 melp_lpc_analysis.m 的输出 a_i 定义，此处需补负号还原
    a_z = [1, a_i]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%matlab自己内置的已经是修改好了的系数，已经是标准化的系数了，所以这里不需要加负号%%%%%%%%%%%
    
    % 仅在调试帧 (41-45帧) 开启详细打印
    is_debug = (frame_counter >= 41 && frame_counter <= 45);
    
    try
        % 3. LPC 转换为弧度 (Radiant) LSF
        lsf_rad = poly2lsf(a_z); 
        
        % 4. 单位转换：Radiant -> Hz -> Scaled Value
        % 转换公式：Hz = Rad * (FS / 2pi); Scaled = Hz * 8.192
        lsf_hz = (lsf_rad' * (FS / (2 * pi)));
        lsf_scaled = lsf_hz * SCALING_FACTOR;
        
        if is_debug
            fprintf('\n[帧%d] LSF特征提取调试:\n', frame_counter);
            fprintf('  原始频率(Hz): [%.1f, %.1f, ..., %.1f]\n', lsf_hz(1), lsf_hz(2), lsf_hz(end));
            fprintf('  码本量级(Scaled): [%.0f, %.0f, ..., %.0f]\n', lsf_scaled(1), lsf_scaled(2), lsf_scaled(end));
        end
        
        % 5. 稳定性与单调性检查 (必须严格递增)
        diff_vals = diff(lsf_scaled);
        if any(diff_vals <= 50) % 间隔需大于一定阈值以防合成滤波器溢出
            if is_debug
                fprintf('  ⚠️ 警告：检测到间隔过小或非单调，触发回退机制。\n');
            end
            lsf = fallback_lsf;
        else
            lsf = lsf_scaled(:); 
        end
        
    catch ME
        if is_debug
            fprintf('  ❌ 转换失败 (%s)，使用回退矢量。\n', ME.message);
        end
        lsf = fallback_lsf;
    end
end