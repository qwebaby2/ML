clear
clc
close
%% OMP主程序
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
At = zeros(M,K);
lambda = zeros(1,K);                
theta_liklihood = zeros(N,1);
for i = 1:K
    inner_product_max = 0;
    index = -1;
    for j = 1:N
        inner_product = abs(A(:,j)' * r);
        if inner_product > inner_product_max
            inner_product_max = inner_product;
            index = j;
        end
    end
    At(:,i) = A(:,index);
    lambda(i) = index;
    theta = ((At(:,1:i)' * At(:,1:i)))^(-1)*At(:,1:i)' * y;
    r = y - At(:,1:i)*theta;
end
for ii = 1:K
    theta_liklihood(lambda(ii)) = theta(ii);
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
