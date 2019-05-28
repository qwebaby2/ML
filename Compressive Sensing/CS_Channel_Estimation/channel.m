function[h]=channel(signal_length,taps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------信道脉冲响应,得h-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taumax=signal_length;
tau=taumax*rand(1,taps); %随机产生延时
tau=ceil(tau);
%path_power=2; %路径功率
ampli_real=randn(1,taps);%根据路径功率用高斯过程得到复抽头系数实部
ampli_img=randn(1,taps);%根据路径功率用高斯过程得到复抽头系数 虚部
ampli=ampli_real+ampli_img*j;
h=zeros(1,taumax);
for i=1:taps
    h(tau(i))=exp(-tau(i)/taumax)*ampli(i)/(abs(ampli(i)));
end
h

