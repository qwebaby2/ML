
%频域信道估计
clc;
clear all;
L1=31;
taps=6;%抽头数
K=taps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------频域的信道脉冲响应----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1=channel(L1,taps);
l2=channel(L1,taps);
h=cat(2,l1,l2)';
L=size(h,1);

cs=zeros(3,7);
ls=zeros(3,7);
mmse=zeros(3,7);

for t=1:3
    N1=16*t;%训练序列长度
N=N1*2;
[cs_mse_ave,ls_mse_ave,mmse_mse_ave]=MSE_com(N,L,K,h,N1);
cs(t,:)=cs_mse_ave;
ls(t,:)=ls_mse_ave;
mmse(t,:)=mmse_mse_ave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------图示比较-----------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%..
for n=0:6
    SNR(n+1)=5*n;
end

figure(1);
semilogy(SNR,cs(1,:),'ro-.','LineWidth',1.5);
grid on;
hold on;
semilogy(SNR,cs(2,:),'rp-','LineWidth',1.5);
hold on;
semilogy(SNR,cs(3,:),'rs-','LineWidth',2.5);
hold on;
semilogy(SNR,ls(2,:),'bo-.','LineWidth',1.5);
hold on;
semilogy(SNR,ls(1,:),'bp-','LineWidth',1.5);
hold on;
semilogy(SNR,ls(3,:),'bs-','LineWidth',2.5);
hold on;
semilogy(SNR,mmse(1,:),'mo-.','LineWidth',1.5);
hold on;
semilogy(SNR,mmse(2,:),'mp-','LineWidth',1.5);
hold on;
semilogy(SNR,mmse(3,:),'ms-','LineWidth',2.5);
hold on;
xlabel('SNR in DB');
ylabel('mean squared error');