clear
clc
close
%% CoSaMP主程序
M = 64;                             % 观测个数
N = 256;                            % 信号长度
phi = randn(M, N);                  % 测量矩阵
psi = eye(N);                       % 稀疏矩阵
A = phi * psi;                      % 传感矩阵
K = 10;                             % 信号稀疏度
Index_K = randperm(N);
x = zeros(N, 1);                    % 传输信号
x(Index_K(1:K)) = 5 * randn(K,1);
y = phi * x;                        % 观测信号
r = y;                              % 初始化残差             
theta_liklihood = zeros(N,1);       % 恢复出的似然投影
inner_product = zeros(N,1);
At = [];
lambda = [];
flag = 1;                           % 迭代标志
iter_times = 0;                     % 迭代次数
epsilon = 1e-7;                     % 迭代允许误差
At_index = [];
while flag
    inner_product = abs(A' * r);
    [inner_product_sorted, inner_product_index] =...
        sort(inner_product, 'descend');
    At_index = union(At_index, inner_product_index(1:2*K));
    At = A(:,At_index); 
    theta = ((At' * At))^(-1) * At' * y;
    [~, theta_index] =sort(abs(theta), 'descend');
    theta_index = sort(theta_index(1:K));
    At = At(:,theta_index);
    At_index = At_index(theta_index);
    r = y - At * theta(theta_index);
    if (r' * r) / M < epsilon
        flag = 0;
    end
    % 迭代停止条件
    iter_times = iter_times + 1;
end
for ii = 1:K
    theta_liklihood(At_index(ii)) = theta(theta_index(ii));
end
x_r = psi * theta_liklihood;
%% 绘图
figure;
plot(x_r,'k.-');                    % 绘出x的恢复信号
hold on;
plot(x,'r');                        % 绘出原信号x
hold off;
legend('Recovery','Original')
fprintf('\n恢复残差：');
norm(x_r-x)                         % 恢复残差
