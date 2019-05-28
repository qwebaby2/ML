clear 
clc       
M = 64;%�۲�ֵ����        
N = 256;%�ź�x�ĳ���        
K = 10;%�ź�x��ϡ���        
Index_K = randperm(N);        
x = zeros(N,1);        
x(Index_K(1:K)) = 5*randn(K,1);%xΪKϡ��ģ���λ��������� 
%x(Index_K(1:K)) = sign(5*randn(K,1));             
Phi = randn(M,N);%��������Ϊ��˹����    
Phi = orth(Phi')';            
sigma = 0.005;      
e = sigma*randn(M,1);    
y = Phi * x + e;%�õ��۲�����y        
% y = Phi * x;%�õ��۲�����y      
%% �ָ��ع��ź�x        
tic
% lamda = sigma*sqrt(2*log(N));
lamda = 0.1*max(abs(Phi'*y));
fprintf('\nlamda = %f\n',lamda)
%x_r =  BPDN_quadprog(y,Phi,lamda); 
x_r = IST_Basic(y,Phi,lamda);          
toc  
%Debias
[xsorted,inds] = sort(abs(x_r), 'descend'); 
AI = Phi(:,inds(xsorted(:)>1e-3));
xI = pinv(AI'*AI)*AI'*y;
x_bias = zeros(length(x),1);
x_bias(inds(xsorted(:)>1e-3)) = xI;
%% ��ͼ        
figure;        
plot(x_r,'k.-');% ���x�Ļָ��ź�        
hold on;        
plot(x,'r');% ���ԭ�ź�x
title("not debias")
hold off;        
legend('Recovery','Original')        
fprintf('\n�ָ��в�(original)��');        
fprintf('%f\n',norm(x_r-x));% �ָ��в�  
% Debias��ƫ��
figure;        
plot(x_bias,'k.-');% ���x�Ļָ��ź�        
hold on;        
plot(x,'r');% ���ԭ�ź�x
title("debias")
hold off;        
legend('Recovery-debise','Original')        
fprintf('�ָ��в�(debias)��');        
fprintf('%f\n',norm(x_bias-x));% �ָ��в�  