clc
clear

load('K_means');
data = X;
K = 3;
clear X;
plot(data(:,1),data(:,2),'k.')
%%%%%%%%%%%%%%%%
[m,n]=size(data);           %��data�Ĵ�С��m���������ĸ�����n��������������
centre_pt=zeros(K,n);       %�洢���ĵ�
distance=zeros(m,K);        %�洢�������������ĵ�ľ��룬�������е�����������ÿ�����������ĵ�ľ���
label=zeros(m,1);           %�洢���������ı�ǩ
iteration_max=500;


for i=1:K
    centre_pt(i,:)=data(floor(rand*m),:);%�������K����������
end



%%���濪ʼ���е���
for iteration=1:iteration_max
    distance1=zeros(m,K);
    for i=1:m
        for j=1:K
            distance1(i,j)=sqrt(sum((data(i,:)-centre_pt(j,:)).^2));%һ�����������о������ĵľ���
        end%һ�д���һ��������K�д�����K���������ĵľ���
    end
[~,label] = min(distance1,[],2);
for i = 1:K
temp_label = find(label == i);
centre_pt(i,:) = sum(data(temp_label,:)) / length(temp_label);
end

myerror=0;
for i=1:m
    for k=1:K%������벻�ٱ仯�������ĵ㲻�ٱ仯�����ߴﵽ�����ĵ�����������ֹͣ����
    myerror=myerror+sum((distance1(i,k)-distance(i,k)).^2);
    end
end
distance=distance1;
%������е����ĵ㲻���ƶ�,����ѭ��
if myerror==0
    break;
end

end
fprintf("һ��������%d��",iteration)
label1 = find(label == 1);
label2 = find(label == 2);
label3 = find(label == 3);
plot(data(label1,1),data(label1,2),'k.');
hold on;
plot(data(label2,1),data(label2,2),'r.');
plot(data(label3,1),data(label3,2),'b.');
hold off;