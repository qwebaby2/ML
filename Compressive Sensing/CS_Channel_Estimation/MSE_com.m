function [cs_mse_ave,ls_mse_ave,mmse_mse_ave]=MSE_com(N,L,K,h,N1)

W_h=1/sqrt(N)*fft(eye(N,L));
H=W_h*h;
H1=H(1:N1,:);
H2=H((N1+1):N,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------training sequence----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=randn(1,N);
d=d/std(d);
d=d-mean(d);
X=diag(d);
X1=X(1:N1,1:N1);
X2=X((N1+1):N,(N1+1):N);
XH=X*H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------求h的自协方差矩阵-Rhh-------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gg=diag(h);
gg_myu = sum(gg, 1)/L;                    
gg_mid = gg - gg_myu(ones(L,1),:);        
sum_gg_mid= sum(gg_mid, 1);
Rgg = (gg_mid' * gg_mid- (sum_gg_mid'  * sum_gg_mid) / L) / (L- 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------添加高斯白噪声，得Y-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1=ones(N,1); 
for m=1:20%多组实验取平均
    for n=0:6
       
SNR(n+1)=5*n;%比较不同SNR
clear j;
n1=n1*0.01j;%保证下面的awgn函数输入的是复高斯噪声
No=awgn(n1,SNR(n+1));%white Gaussian noise
%variance=var(noise);
SNR_log=10^(SNR(n+1)/10);
variance=var(XH)/SNR_log;
No=variance/var(No)*No;
var_No=var(No);
%No=fft(noise);
%Y = AWGN(X,SNR) adds  to X.  The SNR is in dB.The power of X is assumed to be 0 dBW.  If X is complex, then AWGN adds complex noise.
%No=fft(noise);
Y=XH+No;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------LS/MMSE信道估计，得MSE-------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_squared_error_ls=LS_MSE_calc(X,H,Y,N); 
%Evaluating the mean squared error for the MMSE estimator..
mean_squared_error_mmse=MMSE_MSE_calc(X,H,Y,Rgg,var_No,N,L); 
mmse_mse(m,n+1)=mean_squared_error_mmse;
ls_mse(m,n+1)=mean_squared_error_ls;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------CS信道估计H，得MSE--------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS evaluate H
s=Y;
Phi=X;
T=Phi*W_h;                               %  恢复矩阵(测量矩阵*正交反变换矩阵
re_H=zeros(1,N);                         %  待重构的谱域(变换域)向量   
re_y=zeros(1,L);
[pos_arry,aug_y]=omp(K,s,T);              %   pos_arry:最大投影系数对应的位置,
[cos_pos_arry,aug_y]=omp(K,s,T);          %   pos_arry:最大投影系数对应的位置,
re_y(pos_arry)=aug_y;
re_H=W_h*re_y.';                     %  做傅里叶变换重构得到原信号                               

diff_value=abs((re_H) -(H));
re_error=mean((diff_value./abs(H)).^2);
cs_mse(m,n+1)=re_error;
end
end

mmse_mse_ave=mean(mmse_mse);
ls_mse_ave=mean(ls_mse);
cs_mse_ave=mean(cs_mse);