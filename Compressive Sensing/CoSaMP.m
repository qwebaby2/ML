clear
clc
close
%% CoSaMP������
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
theta_liklihood = zeros(N,1);       % �ָ�������ȻͶӰ
inner_product = zeros(N,1);
At = [];
lambda = [];
flag = 1;                           % ������־
iter_times = 0;                     % ��������
epsilon = 1e-7;                     % �����������
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
    % ����ֹͣ����
    iter_times = iter_times + 1;
end
for ii = 1:K
    theta_liklihood(At_index(ii)) = theta(theta_index(ii));
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
