clc
clear

load('K_means');
data = X;
K = 3;
clear X;
plot(data(:,1),data(:,2),'k.')
%%%%%%%%%%%%%%%%
[m,n]=size(data);           %求data的大小，m代表样本的个数，n代表数据特征数
centre_pt=zeros(K,n);       %存储中心点
distance=zeros(m,K);        %存储各个样本到中心点的距离，行是所有的样本，列是每个样本与中心点的距离
label=zeros(m,1);           %存储各个样本的标签
iteration_max=500;


for i=1:K
    centre_pt(i,:)=data(floor(rand*m),:);%随机产生K个聚类中心
end



%%下面开始进行迭代
for iteration=1:iteration_max
    distance1=zeros(m,K);
    for i=1:m
        for j=1:K
            distance1(i,j)=sqrt(sum((data(i,:)-centre_pt(j,:)).^2));%一个样本与所有聚类中心的距离
        end%一行代表一个样本，K列代表与K个聚类中心的距离
    end
[~,label] = min(distance1,[],2);
for i = 1:K
temp_label = find(label == i);
centre_pt(i,:) = sum(data(temp_label,:)) / length(temp_label);
end

myerror=0;
for i=1:m
    for k=1:K%如果距离不再变化，即中心点不再变化，或者达到了最大的迭代次数，则停止迭代
    myerror=myerror+sum((distance1(i,k)-distance(i,k)).^2);
    end
end
distance=distance1;
%如果所有的中心点不再移动,跳出循环
if myerror==0
    break;
end

end
fprintf("一共迭代了%d次",iteration)
label1 = find(label == 1);
label2 = find(label == 2);
label3 = find(label == 3);
plot(data(label1,1),data(label1,2),'k.');
hold on;
plot(data(label2,1),data(label2,2),'r.');
plot(data(label3,1),data(label3,2),'b.');
hold off;