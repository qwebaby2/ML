clear
clc
close
%% OMP������
M = 64;                             % �۲����
N = 256;                            % �źų���
phi = randn(M, N);                  % ��������
psi = eye(N);                       % ϡ�����
A = phi * psi;                      % ���о���
K = 10;                             % �ź�ϡ���
Index_K = randperm(N);
x = zeros(N, 1);                    % �����ź�
x(Index_K(1:K)) = 5 * randn(K,1);
y = phi * x;                        % �۲��ź�
r = y;                              % ��ʼ���в�
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
%% ��ͼ
figure;
plot(x_r,'k.-');                    % ���x�Ļָ��ź�
hold on;
plot(x,'r');                        % ���ԭ�ź�x
hold off;
legend('Recovery','Original')
fprintf('\n�ָ��в');
norm(x_r-x)                         % �ָ��в�
