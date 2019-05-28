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
% --------------------��h����Э�������-Rhh-------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gg=diag(h);
gg_myu = sum(gg, 1)/L;                    
gg_mid = gg - gg_myu(ones(L,1),:);        
sum_gg_mid= sum(gg_mid, 1);
Rgg = (gg_mid' * gg_mid- (sum_gg_mid'  * sum_gg_mid) / L) / (L- 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------��Ӹ�˹����������Y-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1=ones(N,1); 
for m=1:20%����ʵ��ȡƽ��
    for n=0:6
       
SNR(n+1)=5*n;%�Ƚϲ�ͬSNR
clear j;
n1=n1*0.01j;%��֤�����awgn����������Ǹ���˹����
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
%-----------------------LS/MMSE�ŵ����ƣ���MSE-------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_squared_error_ls=LS_MSE_calc(X,H,Y,N); 
%Evaluating the mean squared error for the MMSE estimator..
mean_squared_error_mmse=MMSE_MSE_calc(X,H,Y,Rgg,var_No,N,L); 
mmse_mse(m,n+1)=mean_squared_error_mmse;
ls_mse(m,n+1)=mean_squared_error_ls;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------CS�ŵ�����H����MSE--------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS evaluate H
s=Y;
Phi=X;
T=Phi*W_h;                               %  �ָ�����(��������*�������任����
re_H=zeros(1,N);                         %  ���ع�������(�任��)����   
re_y=zeros(1,L);
[pos_arry,aug_y]=omp(K,s,T);              %   pos_arry:���ͶӰϵ����Ӧ��λ��,
[cos_pos_arry,aug_y]=omp(K,s,T);          %   pos_arry:���ͶӰϵ����Ӧ��λ��,
re_y(pos_arry)=aug_y;
re_H=W_h*re_y.';                     %  ������Ҷ�任�ع��õ�ԭ�ź�                               

diff_value=abs((re_H) -(H));
re_error=mean((diff_value./abs(H)).^2);
cs_mse(m,n+1)=re_error;
end
end

mmse_mse_ave=mean(mmse_mse);
ls_mse_ave=mean(ls_mse);
cs_mse_ave=mean(cs_mse);