% =========================================================================
% 脚本: visualize_encoder_output.m
% 功能: 加载 1200bps 编码器输出 (.mat)，解析 81-bit 码流并可视化核心参数
% =========================================================================
clc; clear; close all;

% --- 1. 加载数据 ---
filename = 'encoder_output_new.mat';
if ~isfile(filename)
    error('找不到文件 %s，请先运行编码器！', filename);
end
fprintf('正在加载 %s ... ', filename);
load(filename); % 预期包含 'all_bit_streams'
fprintf('完成。\n');

if ~exist('all_bit_streams', 'var')
    error('MAT文件中缺少 all_bit_streams 变量！');
end

% --- 2. 解析码流 ---
num_sf = length(all_bit_streams);
fprintf('共发现 %d 个超帧。\n正在解析比特流...\n', num_sf);

% 初始化数据容器
data_gain = zeros(1, num_sf);
data_pitch = zeros(1, num_sf);
data_uv = zeros(1, num_sf);
data_bpvc = zeros(1, num_sf);
data_mode = zeros(1, num_sf); % 0=Mode0, 1=Mode1
data_lsf_main = zeros(1, num_sf); % Mode1:Base / Mode0:F3_S1
data_lsf_sub  = zeros(1, num_sf); % Mode1:Res  / Mode0:F3_S2

for i = 1:num_sf
    b_str = all_bit_streams{i};
    if length(b_str) < 81
        warning('第 %d 帧比特流长度不足 (%d)，跳过。', i, length(b_str));
        continue;
    end
    
    % --- 字段提取 (1-based index) ---
    % 结构: S(1)|UV(3)|Par(1)|Pitch(9)|LSF(42)|Gain(10)|BPVC(6)|FS(8)|J(1)
    
    % 1. UV (Bits 2-4)
    uv_bin = b_str(2:4);
    data_uv(i) = bin2dec(uv_bin);
    
    % 2. Pitch (Bits 6-14)
    pitch_bin = b_str(6:14);
    data_pitch(i) = bin2dec(pitch_bin);
    
    % 3. LSF (Bits 15-56) [42 bits]
    lsf_bin = b_str(15:56);
    
    % 判断模式: UV='001' (Decimal 1) -> Mode 1 (插值)
    if strcmp(uv_bin, '001')
        data_mode(i) = 1; 
        % Mode 1 结构: Base(9) + Int(4) + S1(8) + ...
        % 我们提取 Base 作为主索引，S1 (残差第一级) 作为副索引
        data_lsf_main(i) = bin2dec(lsf_bin(1:9));   % Base Index (0-511)
        data_lsf_sub(i)  = bin2dec(lsf_bin(14:21)); % MSVQ Stage 1 (0-255)
    else
        data_mode(i) = 0;
        % Mode 0 结构: F1(8+6) + F2(8+6) + F3(8+6)
        % 为了对比，我们取 F3 (本超帧最后一帧) 的 S1 作为主索引
        % F3 起始位置: 14 + 14 + 1 = 29
        data_lsf_main(i) = bin2dec(lsf_bin(29:36)); % F3 S1 (0-255)
        data_lsf_sub(i)  = bin2dec(lsf_bin(37:42)); % F3 S2 (0-63)
    end
    
    % 4. Gain (Bits 57-66)
    data_gain(i) = bin2dec(b_str(57:66));
    
    % 5. BPVC (Bits 67-72)
    data_bpvc(i) = bin2dec(b_str(67:72));
end

fprintf('解析完成。Mode 1 帧数: %d / %d\n', sum(data_mode), num_sf);

% --- 3. 可视化绘图 ---
figure('Name', 'MELPe 1200bps Encoder Analysis', 'Position', [100 100 1000 800]);

t_axis = 1:num_sf;

% Subplot 1: Gain
subplot(4,1,1);
plot(t_axis, data_gain, 'b.-', 'LineWidth', 1);
title('① Gain Index (Energy Profile)');
ylabel('Index (0-1023)'); grid on; xlim([1 num_sf]);

% Subplot 2: Pitch
subplot(4,1,2);
plot(t_axis, data_pitch, 'g.-', 'LineWidth', 1);
title('② Pitch Index');
ylabel('Index (0-511)'); grid on; xlim([1 num_sf]);

% Subplot 3: BPVC & Mode
subplot(4,1,3);
bar(t_axis, data_bpvc, 'FaceColor', 'c'); hold on;
% 用红色阶梯线叠加显示 Mode 1 状态
stairs(t_axis, data_mode * 50, 'r-', 'LineWidth', 1.5); 
legend('BPVC Index', 'Mode 1 Flag (Scaled)', 'Location', 'Best');
title('③ BPVC (Cyan Bars) & Coding Mode (Red Line)');
ylabel('Index'); grid on; xlim([1 num_sf]);

% Subplot 4: LSF Analysis
subplot(4,1,4);
hold on;
% 绘制 Mode 0 (普通) 为橙色空心圆
idx_m0 = find(data_mode == 0);
scatter(t_axis(idx_m0), data_lsf_main(idx_m0), 40, 'MarkerEdgeColor', [0.85 0.33 0.1], 'LineWidth', 1.5);
% 绘制 Mode 1 (插值) 为紫色实心点
idx_m1 = find(data_mode == 1);
scatter(t_axis(idx_m1), data_lsf_main(idx_m1), 40, 'MarkerFaceColor', [0.5 0 0.5], 'MarkerEdgeColor', 'none');

plot(t_axis, data_lsf_main, 'k:', 'LineWidth', 0.5); % 连线观察连续性
title('④ LSF Primary Index (Orange=Mode0/F3, Purple=Mode1/Base)');
ylabel('Index'); xlabel('Superframe Index'); 
grid on; xlim([1 num_sf]);
legend('Mode 0 (Independent)', 'Mode 1 (Interpolated)', 'Location', 'Best');

fprintf('✅ 图表绘制完成。\n');