% =========================================================================
% 脚本: visualize_melp_complete.m
% 功能: MELPe 1200bps 终极可视化 (比特流对比 + 参数拆解 + 波形绘制)
% =========================================================================
clc; clear; close all;

fprintf('=== MELPe 1200bps 深度可视化诊断 ===\n');

% --- 1. 加载数据 ---
try
    enc = load('encoder_output_new.mat'); 
    fprintf('✅ 加载编码器数据 (encoder_output_new.mat)\n');
    
    if isfile('decode_output.mat')
        dec = load('decode_output.mat');
        fprintf('✅ 加载解码器数据 (decode_output.mat)\n');
    else
        error('❌ 找不到 decode_output.mat');
    end
catch
    error('❌ 文件加载失败，请确保您已运行过编码器和解码器。');
end

% 获取数据流
enc_streams = enc.all_bit_streams;
dec_streams = [];
if isfield(dec, 'debug_rx_streams')
    dec_streams = dec.debug_rx_streams;
else
    fprintf('⚠️ 警告: 解码器未保存 debug_rx_streams，仅显示编码器流。\n');
end

% 获取信号
sig_orig = []; if isfield(enc, 'signal_final'), sig_orig = enc.signal_final; end
sig_dec  = []; if isfield(dec, 'decoded_signal'), sig_dec = dec.decoded_signal; end

% --- 2. 文本可视化: 比特流与参数表 ---
fprintf('\n');
fprintf('%s\n', repmat('=', 1, 160));
fprintf('FRAME-BY-FRAME BITSTREAM & PARAMETER ANALYSIS\n');
fprintf('%s\n', repmat('=', 1, 160));
fprintf('%-4s | %-6s | %-9s | %-9s | %-9s | %-38s | %-12s | %s\n', ...
    'SF', 'Status', 'Pitch', 'Gain', 'BPVC/FS', 'LSF Breakdown', 'Stream Head', 'Full 81-bit Stream Check');
fprintf('%s\n', repmat('-', 1, 160));

num_frames = length(enc_streams);
p_trace_enc = []; p_trace_dec = [];
g_trace_enc = []; g_trace_dec = [];

for k = 1:num_frames
    % --- 解析编码器流 ---
    es = enc_streams{k};
    if length(es) ~= 81, continue; end
    
    [e_p] = parse_81bits(es);
    
    % --- 解析解码器流 (如果存在) ---
    ds = ''; d_p = [];
    if ~isempty(dec_streams) && k <= length(dec_streams)
        ds = dec_streams{k};
        if length(ds) == 81
            [d_p] = parse_81bits(ds);
        end
    end
    
    % --- 对比 ---
    status = 'MATCH';
    if ~isempty(ds) && ~strcmp(es, ds)
        status = '❌ DIFF'; 
    elseif isempty(ds)
        status = 'NO DEC';
    end
    
    % --- 格式化打印 ---
    % Pitch
    str_pit = sprintf('%d', e_p.pitch);
    if ~isempty(d_p), str_pit = sprintf('%d/%d', e_p.pitch, d_p.pitch); end
    
    % Gain
    str_gain = sprintf('%d', e_p.gain);
    if ~isempty(d_p), str_gain = sprintf('%d/%d', e_p.gain, d_p.gain); end
    
    % BPVC/FS
    str_misc = sprintf('B:%X F:%d', e_p.bpvc, e_p.fs);
    
    % LSF
    if e_p.mode == 1
        str_lsf = sprintf('M1 B:%d I:%d R:%d/%d/%d/%d', e_p.lsf.base, e_p.lsf.int, e_p.lsf.s1, e_p.lsf.s2, e_p.lsf.s3, e_p.lsf.s4);
    else
        str_lsf = sprintf('M0 F1:%d/%d F2:%d/%d F3:%d/%d', e_p.lsf.f1s1, e_p.lsf.f1s2, e_p.lsf.f2s1, e_p.lsf.f2s2, e_p.lsf.f3s1, e_p.lsf.f3s2);
    end
    
    % Bitstream Head
    str_head = sprintf('%s...', es(1:10));
    
    fprintf('%03d  | %-6s | %-9s | %-9s | %-9s | %-38s | %-12s | E: %s\n', ...
        k, status, str_pit, str_gain, str_misc, str_lsf, str_head, es);
        
    if strcmp(status, '❌ DIFF')
        fprintf('      |        |           |           |           |                                        |              | D: %s\n', ds);
        % 标记差异位
        diff_mask = repmat(' ', 1, 81);
        diff_mask(es ~= ds) = '^';
        fprintf('      |        |           |           |           |                                        |              |    %s\n', diff_mask);
    end
    
    % 收集数据画图
    p_trace_enc(end+1) = e_p.pitch;
    g_trace_enc(end+1) = e_p.gain;
end
fprintf('%s\n', repmat('-', 1, 160));

% --- 3. 图形可视化 ---
figure('Name', 'MELPe Complete Visualization', 'Position', [100, 100, 1200, 800]);

% A. 波形对比
subplot(3,1,1);
if ~isempty(sig_orig) && ~isempty(sig_dec)
    % 简单的对齐 (假设无延迟或手动对齐)
    L = min(length(sig_orig), length(sig_dec));
    t = (0:L-1)/8000;
    plot(t, sig_orig(1:L), 'b', 'LineWidth', 1); hold on;
    plot(t, sig_dec(1:L), 'r', 'LineWidth', 1);
    legend('Encoder Input', 'Decoder Output');
    title('Waveform Comparison (Time Domain)');
    grid on;
else
    text(0.5, 0.5, 'No Signal Data Available', 'HorizontalAlignment', 'center');
end

% B. Pitch 轨迹
subplot(3,1,2);
plot(p_trace_enc, 'o-', 'LineWidth', 1.5);
title('Pitch Trajectory (Encoder)');
ylabel('Pitch Index'); grid on;
xlim([1, length(p_trace_enc)]);

% C. Gain 轨迹
subplot(3,1,3);
plot(g_trace_enc, 's-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
title('Gain Trajectory (Encoder)');
ylabel('Gain Index'); grid on;
xlim([1, length(g_trace_enc)]);
xlabel('Superframe Index');

fprintf('✅ 可视化完成！请查看控制台输出表和图形窗口。\n');

% =========================================================================
% 辅助函数: 解析 81-bit 结构
% =========================================================================
function p = parse_81bits(s)
    % 结构: Sync(1)|UV(3)|Par(1)|Pitch(9)|LSF(42)|Gain(10)|BPVC(6)|FS(8)|Jit(1)
    p.sync  = s(1);
    p.uv    = s(2:4);
    p.par   = s(5);
    p.pitch = bin2dec(s(6:14));
    p.lsf_bits = s(15:56);
    p.gain  = bin2dec(s(57:66));
    p.bpvc  = bin2dec(s(67:72));
    p.fs    = bin2dec(s(73:80));
    p.jit   = s(81);
    
    if strcmp(p.uv, '001')
        p.mode = 1;
        p.lsf.base = bin2dec(p.lsf_bits(1:9));
        p.lsf.int  = bin2dec(p.lsf_bits(10:13));
        p.lsf.s1   = bin2dec(p.lsf_bits(14:21));
        p.lsf.s2   = bin2dec(p.lsf_bits(22:27));
        p.lsf.s3   = bin2dec(p.lsf_bits(28:33));
        p.lsf.s4   = bin2dec(p.lsf_bits(34:39));
    else
        p.mode = 0;
        p.lsf.f1s1 = bin2dec(p.lsf_bits(1:8));
        p.lsf.f1s2 = bin2dec(p.lsf_bits(9:14));
        p.lsf.f2s1 = bin2dec(p.lsf_bits(15:22));
        p.lsf.f2s2 = bin2dec(p.lsf_bits(23:28));
        p.lsf.f3s1 = bin2dec(p.lsf_bits(29:36));
        p.lsf.f3s2 = bin2dec(p.lsf_bits(37:42));
    end
end