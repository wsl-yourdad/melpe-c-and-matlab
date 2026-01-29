function a = lsf_to_lpc(lsf, LPC_ORD)
% LSF_TO_LPC  将线谱对(LSF)频率转换回线性预测系数(LPC)。
%
% 输入:
%   lsf:     一个包含LSF频率的向量(单位: Hz)。
%   LPC_ORD: (可选) LPC的阶数，默认为10。
%
% 输出:
%   a:       一个包含LPC系数的列向量 [a1, a2, ..., a_p]'。
%
% 算法:
%   该函数利用LSF的根与单位圆交替出现的特性，构建两个多项式P(z)和Q(z)，
%   然后通过 A(z) = 0.5 * [P(z) + Q(z)] 来恢复LPC多项式A(z)。
%   这是一种数值上非常稳定的标准算法。

    if nargin < 2
        LPC_ORD = 10;
    end
    
    % 确保输入是列向量
    lsf = lsf(:);

    % --- 1. 将LSF频率(Hz)转换为角频率(radians) ---
    % 假设采样率为8000Hz (标准MELP)
    FS = 8000;
    w = lsf * (2 * pi / FS);

    % --- 2. 分离奇数和偶数索引的LSF ---
    % LSFs被分成两组，用于构建P(z)和Q(z)
    p = w(1:2:LPC_ORD); % 奇数索引: w1, w3, ..., w9
    q = w(2:2:LPC_ORD); % 偶数索引: w2, w4, ..., w10

    % P(z) 和 Q(z) 的阶数
    order_p = length(p);
    order_q = length(q);

    % --- 3. 构建P1(z)和Q1(z) ---
    % P(z) = P1(z) * (1+z^-1)
    % Q(z) = Q1(z) * (1-z^-1)
    % P1(z)和Q1(z)的根是 exp(j*w) 和 exp(-j*w)
    
    % 初始化P1和Q1为1 (即 z^0)
    P1 = 1;
    Q1 = 1;

    for i = 1:order_p
        % 每次乘以一个二阶节 (1 - 2*cos(w)*z^-1 + z^-2)
        % 这等同于与 [1, -2*cos(p(i)), 1] 进行卷积
        root_p = [1, -2 * cos(p(i)), 1];
        P1 = conv(P1, root_p);
    end

    for i = 1:order_q
        root_q = [1, -2 * cos(q(i)), 1];
        Q1 = conv(Q1, root_q);
    end

    % --- 4. 构建P(z)和Q(z) ---
    % P(z) = P1(z) * (1+z^-1) -> P1(z) + z^-1 * P1(z)
    P = P1 + [0, P1(1:end-1)];
    
    % Q(z) = Q1(z) * (1-z^-1) -> Q1(z) - z^-1 * Q1(z)
    Q = Q1 - [0, Q1(1:end-1)];
    
    % --- 5. 计算最终的LPC多项式A(z) ---
    % A(z) = 0.5 * [P(z) + Q(z)]
    A = 0.5 * (P + Q);
    
    % --- 6. 格式化输出 ---
    % 返回 a1, a2, ... a_p，不包括a0=1
    a = A(2:end)'; % 确保是列向量

end
