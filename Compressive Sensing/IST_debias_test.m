clear 
clc       
M = 64;%观测值个数        
N = 256;%信号x的长度        
K = 10;%信号x的稀疏度        
Index_K = randperm(N);        
x = zeros(N,1);        
x(Index_K(1:K)) = 5*randn(K,1);%x为K稀疏的，且位置是随机的 
%x(Index_K(1:K)) = sign(5*randn(K,1));             
Phi = randn(M,N);%测量矩阵为高斯矩阵    
Phi = orth(Phi')';            
sigma = 0.005;      
e = sigma*randn(M,1);    
y = Phi * x + e;%得到观测向量y        
% y = Phi * x;%得到观测向量y      
%% 恢复重构信号x        
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
%% 绘图        
figure;        
plot(x_r,'k.-');% 绘出x的恢复信号        
hold on;        
plot(x,'r');% 绘出原信号x
title("not debias")
hold off;        
legend('Recovery','Original')        
fprintf('\n恢复残差(original)：');        
fprintf('%f\n',norm(x_r-x));% 恢复残差  
% Debias除偏后
figure;        
plot(x_bias,'k.-');% 绘出x的恢复信号        
hold on;        
plot(x,'r');% 绘出原信号x
title("debias")
hold off;        
legend('Recovery-debise','Original')        
fprintf('恢复残差(debias)：');        
fprintf('%f\n',norm(x_bias-x));% 恢复残差  