%% ========================================================================
%  音频预处理工具: 8kHz采样率 + 5秒定长 + 单声道
%  适用场景: 语音编解码测试 (MELP/LPC 等)
% =========================================================================
clc; clear; close all;

% --- 1. 参数设置 ---
TARGET_FS = 8000;      % 目标采样率 8kHz
TARGET_DUR = 5.0;      % 目标时长 5秒
OUTPUT_NAME = 'input_8k_5s.wav'; % 输出文件名

% --- 2. 选择音频文件 ---
fprintf('请选择一个音频文件...\n');
[file, path] = uigetfile({'*.wav;*.mp3;*.m4a;*.flac;*.aac', '音频文件'}, ...
                         '选择要处理的音频');

if isequal(file, 0)
    disp('用户取消了操作');
    return;
end

input_full_path = fullfile(path, file);
fprintf('正在读取: %s\n', file);

% --- 3. 读取音频 ---
[x, fs_orig] = audioread(input_full_path);

% --- 4. 转单声道 (混合左右声道) ---
% 编解码器通常只需要单声道
if size(x, 2) > 1
    x = mean(x, 2); 
    fprintf('检测到立体声，已合并为单声道。\n');
else
    fprintf('输入已是单声道。\n');
end

% --- 5. 重采样 (Resample) 到 8000Hz ---
if fs_orig ~= TARGET_FS
    fprintf('正在重采样: %d Hz -> %d Hz\n', fs_orig, TARGET_FS);
    % 使用 MATLAB 自带的 polyphase 重采样，抗混叠效果好
    x_resampled = resample(x, TARGET_FS, fs_orig);
else
    x_resampled = x;
    fprintf('采样率已匹配，无需重采样。\n');
end

% --- 6. 截断 或 补零 (定长 5秒) ---
target_samples = round(TARGET_FS * TARGET_DUR);
current_samples = length(x_resampled);

if current_samples > target_samples
    % 情况A: 太长了，直接截断
    x_final = x_resampled(1:target_samples);
    fprintf('音频过长，已截取前 5 秒。\n');
    
elseif current_samples < target_samples
    % 情况B: 太短了，补零 (静音)
    padding = zeros(target_samples - current_samples, 1);
    x_final = [x_resampled; padding];
    fprintf('音频不足 5 秒，已末尾补零。\n');
    
else
    % 情况C: 刚好
    x_final = x_resampled;
end

% --- 7. 幅度归一化 (可选，防止爆音) ---
% 缩放到 -0.99 到 0.99 之间
max_val = max(abs(x_final));
if max_val > 0
    x_final = x_final / max_val * 0.99;
end

% --- 8. 保存文件 ---
audiowrite(OUTPUT_NAME, x_final, TARGET_FS);
fprintf('✅ 处理完成！已保存为: %s\n', OUTPUT_NAME);

% --- 9. 画图验证 ---
t = (0:length(x_final)-1) / TARGET_FS;

figure('Name', '音频处理结果验证', 'Color', 'w');
subplot(2,1,1);
plot(x); title(['原始信号 (', num2str(fs_orig), ' Hz)']); grid on;
axis tight;

subplot(2,1,2);
plot(t, x_final, 'r'); title(['处理后 (8kHz, 5s, 单声道)']); grid on;
xlabel('时间 (秒)');
axis([0 5 -1 1]); % 强制显示 0-5秒 范围

% 播放处理后的声音 (取消注释可直接听到)
% sound(x_final, TARGET_FS);