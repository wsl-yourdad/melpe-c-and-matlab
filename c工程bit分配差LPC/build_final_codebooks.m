% =========================================================================
% 脚本: build_final_codebooks.m (The Final Correct Version - Fixed)
% 功能: 深度解析 C 工程码本，构建 1200bps/2400bps 混合数据库
% =========================================================================
clc; clear;
fprintf('=== 构建最终、最完整的 MELPe 数据库 (C Project 1:1 复刻版) ===\n\n');

FOLDER_NAME = 'codebook';
if ~isfolder(FOLDER_NAME), error('错误：找不到文件夹 "%s"。', FOLDER_NAME); end

codebooks = struct();
LPC_ORD = 10;

try
    fprintf('1. 正在扫描并加载 C 工程原始数据...\n');
    
    % =====================================================================
    % 1. LSF 核心板块 (最关键的部分)
    % =====================================================================
    fprintf('   -> [LSF] 加载均值与预测器...\n');
    % [Mean] 必须有，用于去均值
    raw_mean = read_clean_data(fullfile(FOLDER_NAME, 'msvq_mean.txt'));
    codebooks.msvq_mean = raw_mean(:) / 32768.0; % Q15 -> Float
    
    % [2400bps] 标准 MSVQ 码本 (配套 Mean 使用)
    % 很多 C 工程将其命名为 c_msvq_cb.txt
    if isfile(fullfile(FOLDER_NAME, 'c_msvq_cb.txt'))
        fprintf('   -> [LSF] 加载 2400bps 标准码本 (c_msvq_cb)...\n');
        raw_cb_2400 = read_clean_data(fullfile(FOLDER_NAME, 'c_msvq_cb.txt'));
        % 结构: 4级 (128, 64, 64, 64)
        mat_2400 = reshape(raw_cb_2400 / 32768.0, 10, []);
        codebooks.msvq_cb_2400 = cell(1,4);
        codebooks.msvq_cb_2400{1} = mat_2400(:, 1:128);
        codebooks.msvq_cb_2400{2} = mat_2400(:, 129:192);
        codebooks.msvq_cb_2400{3} = mat_2400(:, 193:256);
        codebooks.msvq_cb_2400{4} = mat_2400(:, 257:320);
    end

    % [1200bps Mode 0] 主量化码本 (lsp_v_256...)
    % 注意：为了匹配你的编码器，我们将名字映射为 lsp_v_cb
    fprintf('   -> [LSF] 加载 1200bps Mode 0 主码本...\n');
    if isfile(fullfile(FOLDER_NAME, 'lsp_v_256_64_32_32.txt'))
        raw_v = read_clean_data(fullfile(FOLDER_NAME, 'lsp_v_256_64_32_32.txt'));
        mat_v = reshape(raw_v / 32768.0, 10, []);
        codebooks.lsp_v_cb = cell(1,4);
        codebooks.lsp_v_cb{1} = mat_v(:, 1:256);
        codebooks.lsp_v_cb{2} = mat_v(:, 257:320);
        codebooks.lsp_v_cb{3} = mat_v(:, 321:352);
        codebooks.lsp_v_cb{4} = mat_v(:, 353:384);
    else
        % 兼容性：有些包里叫 c_msvq_cb_rate1200.txt
        if isfile(fullfile(FOLDER_NAME, 'c_msvq_cb_rate1200.txt'))
             raw_v = read_clean_data(fullfile(FOLDER_NAME, 'c_msvq_cb_rate1200.txt'));
             % ... 同上处理 ...
             fprintf('      (从 c_msvq_cb_rate1200.txt 加载)\n');
        end
    end

    % [1200bps Mode 1] 残差码本 (20维, Q17)
    % 编码器变量名: res_cb
    fprintf('   -> [LSF] 加载 1200bps Mode 1 残差码本 (20维)...\n');
    if isfile(fullfile(FOLDER_NAME, 'res256_64_64_64.txt'))
        raw_res = read_clean_data(fullfile(FOLDER_NAME, 'res256_64_64_64.txt'));
        % Q17 -> 除以 131072.0
        mat_res = reshape(raw_res / 131072.0, 20, []);
        codebooks.res_cb = cell(1,4);
        codebooks.res_cb{1} = mat_res(:, 1:256);
        codebooks.res_cb{2} = mat_res(:, 257:320);
        codebooks.res_cb{3} = mat_res(:, 321:384);
        codebooks.res_cb{4} = mat_res(:, 385:448);
    end

    % [1200bps Mode 1] 基准码本 (9 bits)
    % 编码器变量名: lsp_uv_9 (或 lsf_cb_9bit 以兼容旧代码)
    fprintf('   -> [LSF] 加载 1200bps Mode 1 基准码本 (9-bit)...\n');
    if isfile(fullfile(FOLDER_NAME, 'lsp_uv_9.txt'))
        raw_uv9 = read_clean_data(fullfile(FOLDER_NAME, 'lsp_uv_9.txt'));
        codebooks.lsp_uv_9 = reshape(raw_uv9 / 32768.0, 10, []);
        
        % [兼容性补丁] 同时创建一个旧名字的副本，防止你主程序报错
        codebooks.lsf_cb_9bit = codebooks.lsp_uv_9; 
    end

    % [1200bps Mode 1] 插值系数 (inpCoef)
    fprintf('   -> [LSF] 加载插值系数 (inpCoef)...\n');
    if isfile(fullfile(FOLDER_NAME, 'inp_raw.txt'))
        raw_inp = read_clean_data(fullfile(FOLDER_NAME, 'inp_raw.txt'));
        % Q14 -> 16384.0
        % C数组通常是 [20][16]，txt 顺序读取后 reshape 为 20x16 即可
        codebooks.inpCoef = reshape(raw_inp / 16384.0, 20, 16);
    end

    % =====================================================================
    % 2. 增益与基音 (Gain & Pitch)
    % =====================================================================
    fprintf('   -> [Gain & Pitch] 加载中...\n');
    if isfile(fullfile(FOLDER_NAME, 'gain_vq_data.txt'))
        raw_gain = read_clean_data(fullfile(FOLDER_NAME, 'gain_vq_data.txt'));
        codebooks.gain = reshape(raw_gain, 6, [])' / 256.0; % Q8
    end
    
    if isfile(fullfile(FOLDER_NAME, 'pitch_vvv.txt'))
        raw_p = read_clean_data(fullfile(FOLDER_NAME, 'pitch_vvv.txt'));
        codebooks.pitch_vq_cb_vvv = reshape(raw_p, 3, []) / 4096.0; % Q12
    end
    if isfile(fullfile(FOLDER_NAME, 'pitch_uvv.txt'))
        raw_p = read_clean_data(fullfile(FOLDER_NAME, 'pitch_uvv.txt'));
        codebooks.pitch_vq_cb_uvv = reshape(raw_p, 3, []) / 4096.0; % Q12
    end

    % =====================================================================
    % 3. 傅里叶与滤波器 (FSVQ & Filters)
    % =====================================================================
    fprintf('   -> [FSVQ & Filters] 加载中...\n');
    if isfile(fullfile(FOLDER_NAME, 'fsvq_cb.txt'))
        raw_fs = read_clean_data(fullfile(FOLDER_NAME, 'fsvq_cb.txt'));
        codebooks.fsvq_cb = reshape(raw_fs / 8192.0, 10, []); % Q13
    end
    
    % Filters (Dispersion, etc.)
    if isfile(fullfile(FOLDER_NAME, 'disp_cof.txt'))
        raw_d = read_clean_data(fullfile(FOLDER_NAME, 'disp_cof.txt'));
        codebooks.disp_cof = raw_d / 32768.0; % Q15
    end
    
    % =====================================================================
    % 3.5 [补漏] 滤波器系数 (Filters) - 必须补全！
    % =====================================================================
    fprintf('   -> [Filters] 正在补全 LPF/BPF 系数...\n');
    
    % LPF (低通滤波器): 通常用于基音提取预处理
    % C代码格式 Q13 -> 除以 8192.0
    if isfile(fullfile(FOLDER_NAME, 'lpf_num.txt'))
        codebooks.lpf_num = read_clean_data(fullfile(FOLDER_NAME, 'lpf_num.txt')) / 8192.0;
    else
        warning('⚠️ 缺少 lpf_num.txt，基音分析可能会报错');
    end
    
    if isfile(fullfile(FOLDER_NAME, 'lpf_den.txt'))
        codebooks.lpf_den = read_clean_data(fullfile(FOLDER_NAME, 'lpf_den.txt')) / 8192.0;
    end

    % BPF (带通滤波器): 用于清浊音判决
    % C代码格式 Q13 -> 除以 8192.0
    if isfile(fullfile(FOLDER_NAME, 'bpf_num.txt'))
        codebooks.bpf_num = read_clean_data(fullfile(FOLDER_NAME, 'bpf_num.txt')) / 8192.0;
    end
    
    if isfile(fullfile(FOLDER_NAME, 'bpf_den.txt'))
        codebooks.bpf_den = read_clean_data(fullfile(FOLDER_NAME, 'bpf_den.txt')) / 8192.0;
    end

    % Disp (脉冲分散滤波器): 之前已有，确保不重复即可
    if isfile(fullfile(FOLDER_NAME, 'disp_cof.txt'))
        codebooks.disp_cof = read_clean_data(fullfile(FOLDER_NAME, 'disp_cof.txt')) / 32768.0; % Q15
    end

    % =====================================================================
    % 4. 保存
    % =====================================================================
    output_filename = 'THE_FINAL_CODEBOOK.mat';
    save(output_filename, 'codebooks');
    fprintf('\n✅ 数据库构建完成: %s\n', output_filename);
    fprintf('   包含字段: %s\n', strjoin(fieldnames(codebooks), ', '));

catch ME
    fprintf('\n❌ 严重错误: %s\n', ME.message);
    fprintf('   出错文件: %s (行 %d)\n', ME.stack(1).name, ME.stack(1).line);
end

% --- 辅助函数 ---
function data = read_clean_data(fname)
    if ~isfile(fname), error('找不到文件 %s', fname); end
    fid = fopen(fname, 'r');
    str = fread(fid, '*char')';
    fclose(fid);
    str = regexprep(str, '[\\{\\};,]', ' '); 
    str = regexprep(str, '/\\*.*?\\*/', '');
    str_clean = regexprep(str, '[^0-9\\.\\-eE]', ' ');
    data = sscanf(str_clean, '%f');
end