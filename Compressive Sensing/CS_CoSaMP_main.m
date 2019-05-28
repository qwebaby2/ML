%ѹ����֪�ع��㷨����
clear
close
clc

M = 64;%�۲�ֵ����
N = 256;%�ź�x�ĳ���
K = 12;%�ź�x��ϡ���
Index_K = randperm(N);
x = zeros(N,1);
x(Index_K(1:K)) = 5*randn(K,1);%xΪKϡ��ģ���λ���������
Psi = eye(N);%x������ϡ��ģ�����ϡ�����Ϊ��λ��x=Psi*theta
Phi = randn(M,N);%��������Ϊ��˹����
A = Phi * Psi;%���о���
y = Phi * x;%�õ��۲�����y
%% �ָ��ع��ź�x
tic
theta = CS_CoSaMP( y,A,K );
x_r = Psi * theta;% x=Psi * theta
toc
%% ��ͼ
figure;
plot(x_r,'k.-');%���x�Ļָ��ź�
hold on;
plot(x,'r');%���ԭ�ź�x
hold off;
legend('Recovery','Original')
fprintf('\n�ָ��в');
norm(x_r-x)%�ָ��в�
