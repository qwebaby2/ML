function[h]=channel(signal_length,taps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------�ŵ�������Ӧ,��h-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taumax=signal_length;
tau=taumax*rand(1,taps); %���������ʱ
tau=ceil(tau);
%path_power=2; %·������
ampli_real=randn(1,taps);%����·�������ø�˹���̵õ�����ͷϵ��ʵ��
ampli_img=randn(1,taps);%����·�������ø�˹���̵õ�����ͷϵ�� �鲿
ampli=ampli_real+ampli_img*j;
h=zeros(1,taumax);
for i=1:taps
    h(tau(i))=exp(-tau(i)/taumax)*ampli(i)/(abs(ampli(i)));
end
h

