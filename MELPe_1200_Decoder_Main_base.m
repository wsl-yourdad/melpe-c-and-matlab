%% =========================================================================
% MELPe 1200bps 解码器主程序 (MELPe_1200_Decoder_Main.m)
% 
% 对应 C工程: melpe_s() 和 main() 的逻辑
% =========================================================================
clc; clear; close all;

%% 1. 加载编码器输出和码本
fprintf('=== MELPe 1200bps 解码器启动 ===\n');

% 加载编码后的比特流 (变量名: all_bit_streams)
if ~isfile('encoder_output.mat')
   error('❌ 找不到 encoder_output.mat，请先运行编码器！');
end
load('encoder_output.mat'); 
fprintf('✅ 已加载比特流，共 %d 个超帧\n', length(all_bit_streams));

% 加载基础码本
codebook_file = 'THE_FINAL_CODEBOOK.mat';
if ~isfile(codebook_file)
    error('❌ 找不到码本文件 %s！', codebook_file);
end
load(codebook_file); % 载入 codebooks 结构体

% ...
load('THE_FINAL_CODEBOOK.mat'); 
% [关键] 检查码本维度转置
if isfield(codebooks, 'lsf_cb_mode_C')
    [r, c] = size(codebooks.lsf_cb_mode_C{1});
    if r > c && c == 10
        fprintf('⚠ 检测到 LSP 码本转置，正在自动修正...\n');
        for k = 1:4, codebooks.lsf_cb_mode_C{k} = codebooks.lsf_cb_mode_C{k}'; end
    end
    codebooks.lsp_vq = codebooks.lsf_cb_mode_C;
end
% ...

% =========================================================================
% [关键修正] 码本字段名映射 (Mapping)
% 你的码本生成脚本把 MSVQ 码本叫 'lsf_cb_mode_C'，解码器需要 'lsp_vq'
% =========================================================================
if isfield(codebooks, 'lsf_cb_mode_C')
    codebooks.lsp_vq = codebooks.lsf_cb_mode_C;
    fprintf('✅ 已自动映射: lsf_cb_mode_C -> lsp_vq\n');
elseif isfield(codebooks, 'lsf_cb')
    codebooks.lsp_vq = codebooks.lsf_cb;
    fprintf('✅ 已自动映射: lsf_cb -> lsp_vq\n');
else
    % 如果连 mode_C 都没有，可能需要检查变量名
    disp('⚠ 警告: 未找到 lsf_cb_mode_C，尝试列出所有字段:');
    disp(fieldnames(codebooks));
end

%% ... (前面加载 encoder_output.mat 的代码不变) ...

% 加载码本
load('THE_FINAL_CODEBOOK.mat'); 

% =========================================================================
% [致命错误修复] 强制检查码本维度 (Transposition Check)
% 原因: 如果码本是 Nx10 而不是 10xN，取出来的 LSF 全是错的，导致针状波形。
% =========================================================================
fprintf('正在检查码本维度...\n');

% 1. 检查 lsf_cb_mode_C (LSP VQ)
if isfield(codebooks, 'lsf_cb_mode_C')
    % 检查第一级码本
    [r, c] = size(codebooks.lsf_cb_mode_C{1});
    if r > c && c == 10
        fprintf('⚠ 检测到 LSP 码本转置 (Size: %dx%d)，正在修正...\n', r, c);
        for k = 1:4
            codebooks.lsf_cb_mode_C{k} = codebooks.lsf_cb_mode_C{k}';
        end
    end
    % 映射标准名称
    codebooks.lsp_vq = codebooks.lsf_cb_mode_C;
end

% 2. 检查 9-bit 码本 (Mode 1 用)
if isfield(codebooks, 'lsf_cb_9bit')
    [r, c] = size(codebooks.lsf_cb_9bit);
    if r > c && c == 10
        fprintf('⚠ 检测到 9-bit LSF 码本转置，正在修正...\n');
        codebooks.lsf_cb_9bit = codebooks.lsf_cb_9bit';
    end
end

% 3. 检查 FSVQ 码本
if isfield(codebooks, 'fsvq_vq')
    [r, c] = size(codebooks.fsvq_vq);
    if r > c && c == 10
        fprintf('⚠ 检测到 FSVQ 码本转置，正在修正...\n');
        codebooks.fsvq_vq = codebooks.fsvq_vq';
    end
    codebooks.fsvq_cb = codebooks.fsvq_vq;
end

% 4. 检查 Pitch VQ
if isfield(codebooks, 'pitch_vq_cb_vvv')
    [r, c] = size(codebooks.pitch_vq_cb_vvv);
    if r > c && c == 3
        fprintf('⚠ 检测到 Pitch VQ 码本转置，正在修正...\n');
        codebooks.pitch_vq_cb_vvv = codebooks.pitch_vq_cb_vvv';
    end
    codebooks.pitch_vvv = codebooks.pitch_vq_cb_vvv;
end

% ... (后续的手动注入 res_cb 代码保持不变) ...

% =========================================================================
% [关键修正] 手动注入 20维残差码本 (res_cb)
% 逻辑源自你的 Word 文档 "编码器主循环" 部分
% =========================================================================
if ~isfield(codebooks, 'res_cb') || isempty(codebooks.res_cb)
    fprintf('⚠ 正在手动注入残差码本 (res256_64_64_64.txt)...\n');
    
    res_txt_file = 'res256_64_64_64.txt';
    if ~isfile(res_txt_file)
        error('❌ 致命错误: 找不到 %s 文件！请确保它和脚本在同一目录下。', res_txt_file);
    end
    
    fid = fopen(res_txt_file, 'r');
    % 稳健读取：读成字符串，去除非数字字符，再转数字
    content = fread(fid, '*char')'; 
    fclose(fid);
    
    % 清洗数据：把逗号、斜杠等非数字符号变为空格
    content(content==',' | content=='/' | content=='*' | content=='Q') = ' ';
    raw_data = sscanf(content, '%f');
    
    % Q17 -> Q15 并重塑为 20 维 (448个向量)
    raw_data = raw_data / 4.0;
    
    % 检查数据长度是否足够
    if length(raw_data) < 8960
         % 有时候最后会有多余字符，取前 8960 个即可 (20 * 448)
         if length(raw_data) < 20*448
             error('❌ 残差码本数据长度不足！预期 8960，实际 %d', length(raw_data));
         end
    end
    
    RES_FULL = reshape(raw_data(1:8960), 20, 448);
    
    % 切分并存入 codebooks.res_cb 字段
    codebooks.res_cb = cell(1,4);
    codebooks.res_cb{1} = RES_FULL(:, 1:256);   % Stage 1 (7 bits)
    codebooks.res_cb{2} = RES_FULL(:, 257:320); % Stage 2 (6 bits)
    codebooks.res_cb{3} = RES_FULL(:, 321:384); % Stage 3 (6 bits)
    codebooks.res_cb{4} = RES_FULL(:, 385:448); % Stage 4 (6 bits)
    
    fprintf('✅ res_cb 字段已成功注入 (20x448)。\n');
end

%% 2. 初始化解码器状态 (State Initialization)
% 初始化上一帧参数 (初始为静音/背景噪)
prev_params.lsf = [0.1:0.04:0.46]' * pi; % 初始LSF分布
prev_params.pitch = 50;                  % 初始基音
prev_params.log_pitch = log10(50);
prev_params.gain = [10; 10];             % 初始增益(很小)
prev_params.uv_flag = 1;                 % 初始为浊音
prev_params.jitter = 0;

% 初始化滤波器状态 (LPC合成滤波器记忆)
lpc_filter_state = zeros(10, 1);  % 10阶LPC
% 初始化脉冲相位
pulse_phase = 0;

% 输出音频缓存
final_speech = [];
FRAME_LEN = 180; % 22.5ms @ 8000Hz

%% 3. 解码主循环 (Processing Loop)
fprintf('正在解码...\n');

for s_idx = 1 : length(all_bit_streams)
    
    % --- A. 获取当前81位比特流 ---
    bit_string = all_bit_streams{s_idx}; 
    if iscell(bit_string), bit_string = bit_string{1}; end
    
    % --- B. 解包与参数恢复 ---
    % 调用我们刚写的 unpack 函数
    [params_batch, prev_params] = melp_unpack_1200(bit_string, prev_params, codebooks);
    
    % --- C. 依次合成3个子帧 ---
    for f = 1 : 3
        curr_frame_params = params_batch(f);
        
        % 调用单帧合成器 (注意传入并更新 pulse_phase)
        [pcm_out, lpc_filter_state, pulse_phase] = melp_syn_1200(curr_frame_params, lpc_filter_state, pulse_phase, codebooks);
        
        % 拼接音频
        final_speech = [final_speech; pcm_out];
    end
    
    if mod(s_idx, 10) == 0
        fprintf('.');
    end
end
fprintf('\n✅ 解码完成！\n');

%% 4. 后处理与播放
% 去掉直流偏置
final_speech = final_speech - mean(final_speech);
% 归一化防爆音
max_val = max(abs(final_speech));
if max_val > 0
    final_speech = final_speech / max_val * 0.9;
end

output_filename = 'decoded_output_1200.wav';
audiowrite(output_filename, final_speech, 8000);
fprintf('✅ 音频已保存: %s\n', output_filename);

% 尝试播放
soundsc(final_speech, 8000);

% 画图对比
figure('Name', 'MELPe 1200 Decoding Result');
subplot(2,1,1); 
plot(final_speech); 
title('Decoded Speech (1200 bps)'); 
grid on;

if isfile('test_5s.wav')
    [orig_sig, fs_orig] = audioread('test_5s.wav');
    % 简单的重采样对齐用于画图（如果采样率不同）
    if fs_orig ~= 8000
        orig_sig = resample(orig_sig, 8000, fs_orig);
    end
    % 截取相同长度
    min_len = min(length(orig_sig), length(final_speech));
    subplot(2,1,2); 
    plot(orig_sig(1:min_len)); 
    title('Original Speech'); 
    grid on;
end