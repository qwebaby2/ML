clear
close
clc
%% Import data
data = csvread('input.csv');
features = data(:,1:-1);
%% Parameters of SVM
alphas = zeros(num_data,1);     % lagrangian multiplier
b = 0;
error = zeros(num_data,2);
tol = 0.001;
C = 0.5;
sig = 18;
iter = 0;
max_iter = 10;
%%
alpha_change = 0;
entireSet = 1;%��Ϊһ����ǿ���ѡ��ȫ�������ǲ��ֱ���
while (iter < max_iter) && ((alpha_change > 0) || entireSet)
    alpha_change = 0;
    %% -----------ȫ��������-------------------------
    if entireSet 
        for i = 1:num_data
            Ei = kernel(data,alphas,label,b,i,sig);%�������
            if (label(i)*Ei<-0.001 && alphas(i)<C)||...
                    (label(i)*Ei>0.001 && alphas(i)>0)
                %ѡ����һ��alphas
                [j,Ej] = select(i,data,num_data,alphas,label,b,C,Ei,entireSet,sig);
                alpha_I_old = alphas(i);
                alpha_J_old = alphas(j);
                if label(i) ~= label(j)
                    L = max(0,alphas(j) - alphas(i));
                    H = min(C,C + alphas(j) - alphas(i));
                else
                    L = max(0,alphas(j) + alphas(i) -C);
                    H = min(C,alphas(j) + alphas(i));
                end
                if L==H
                    continue;
                end
                eta = 2*exp(-1/(2*sig)*norm(data(i,:)-data(j,:))^2)- exp(-1/(2*sig)*norm(data(i,:)-data(i,:))^2)...
                     - exp(-1/(2*sig)*norm(data(j,:)-data(j,:))^2);
                if eta >= 0 
                    continue;
                end
                alphas(j) = alphas(j) - label(j)*(Ei-Ej)/eta;
                %���Ʒ�Χ
                if alphas(j) > H
                    alphas(j) = H;
                elseif alphas(j) < L
                    alphas(j) = L;
                end
                if abs(alphas(j) - alpha_J_old) < 1e-4
                    continue;
                end
                alphas(i) = alphas(i) + label(i)*...
                    label(j)*(alpha_J_old-alphas(j));
                b1 = b - Ei - label(i)*(alphas(i)-alpha_I_old)*...
                    exp(-1/(2*sig)*norm(data(i)-data(i))^2)- label(j)*...
                    (alphas(j)-alpha_J_old)*exp(-1/(2*sig)*norm(data(i)-data(j))^2);
                b2 = b - Ej - label(i)*(alphas(i)-alpha_I_old)*...
                    exp(-1/(2*sig)*norm(data(i)-data(j))^2)- label(j)*...
                    (alphas(j)-alpha_J_old)*exp(-1/(2*sig)*norm(data(j)-data(j))^2);
                if (alphas(i) > 0) && (alphas(i) < C)
                    b = b1;
                elseif (alphas(j) > 0) && (alphas(j) < C)
                    b = b2;
                else
                    b = (b1+b2)/2;
                end
                alpha_change = alpha_change + 1;
            end
        end
         iter = iter + 1;
   %% --------------���ֱ���(alphas=0~C)������--------------------------
    else
        index = find(alphas>0 & alphas < C);
        for ii = 1:length(index)
            i = index(ii);
            Ei = kernel(data,alphas,label,b,i,sig);%�������
            if (label(i)*Ei<-0.001 && alphas(i)<C)||...
                    (label(i)*Ei>0.001 && alphas(i)>0)
                %ѡ����һ������
                [j,Ej] = select(i,data,num_data,alphas,label,b,C,Ei,entireSet,sig);
                alpha_I_old = alphas(i);
                alpha_J_old = alphas(j);
                if label(i) ~= label(j)
                    L = max(0,alphas(j) - alphas(i));
                    H = min(C,C + alphas(j) - alphas(i));
                else
                    L = max(0,alphas(j) + alphas(i) -C);
                    H = min(C,alphas(j) + alphas(i));
                end
                if L==H
                    continue;end
                eta = 2*exp(-1/(2*sig)*norm(data(i,:)-data(j,:))^2)- exp(-1/(2*sig)*norm(data(i,:)-data(i,:))^2)...
                    - exp(-1/(2*sig)*norm(data(j,:)-data(j,:))^2);
                if eta >= 0
                    continue;end
                alphas(j) = alphas(j) - label(j)*(Ei-Ej)/eta;  
                %���Ʒ�Χ
                if alphas(j) > H
                    alphas(j) = H;
                elseif alphas(j) < L
                    alphas(j) = L;
                end
                if abs(alphas(j) - alpha_J_old) < 1e-4
                    continue;end
                alphas(i) = alphas(i) + label(i)*...
                    label(j)*(alpha_J_old-alphas(j));
                b1 = b - Ei - label(i)*(alphas(i)-alpha_I_old)*...
                    exp(-1/(2*sig)*norm(data(i,:)-data(i,:))^2)- label(j)*...
                    (alphas(j)-alpha_J_old)*exp(-1/(2*sig)*norm(data(i,:)-data(j,:))^2);
                b2 = b - Ej - label(i)*(alphas(i)-alpha_I_old)*...
                    exp(-1/(2*sig)*norm(data(i,:)-data(j,:))^2)- label(j)*...
                    (alphas(j)-alpha_J_old)*exp(-1/(2*sig)*norm(data(j,:)-data(j,:))^2);
                if (alphas(i) > 0) && (alphas(i) < C)
                    b = b1;
                elseif (alphas(j) > 0) && (alphas(j) < C)
                    b = b2;
                else
                    b = (b1+b2)/2;
                end
                alpha_change = alpha_change + 1;
            end
        end
        iter = iter + 1;
    end
    %% --------------------------------
    if entireSet %��һ��ȫ�����ˣ���һ�ξͱ�ɲ��ֱ���
        entireSet = 0;
    elseif alpha_change == 0 
        %������ֱ������ж�û���ҵ���Ҫ������alpha���ٸ�Ϊȫ����
        entireSet = 1;
    end
    disp(['iter ================== ',num2str(iter)]);    
end
%% ����ȨֵW
predict = zeros(N,1);
W = (alphas.*label)'*data;
%��¼֧������λ��
index_sup = find(alphas ~= 0);
%����Ԥ����
Kernel = zeros(N,N);
for i_outside =1:N
    for j_inside = 1:N
        Kernel(i_outside,j_inside) =exp(-1/(2*sig)*norm(data(i_outside,:)-data(j_inside,:))^2);
    end
end
for k =1:N
predict(k) = (alphas.*label)'*Kernel(:,k);
end
predict = sign(predict);
%% ��ʾ���
figure;
index1 = find(predict==-1);
data1 = (data(index1,:))';
plot3(data1(1,:),data1(2,:),data1(3,:),'b*');
hold on
index2 = find(predict==1);
data2 = (data(index2,:))';
plot3(data2(1,:),data2(2,:),data2(3,:),'+r');
hold on
dataw = (data(index_sup,:))';
plot3(dataw(1,:),dataw(2,:),dataw(3,:),'go','LineWidth',0.5);
title(['�ɳڱ�����ΧC = ',num2str(C)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������֤��������ģ��׼ȷ�ʡ�
NS = 4000;                                           %��֤����Ŀ
varify_data = zeros(3,NS);
varify_label = zeros(1,NS);
correct = 0;
pre_L = zeros(NS,1);
for i=1:NS                                           %at the beginning of every round, to produce the training set in random;
                                                     %make the rand vectors X;
    if  rand>0.5                                     %make sure the points occur in random. and half in half for each class;
        theta=pi*rand;                               %set the angle;
        lenth=r+wth*(rand-0.5);                      %the range;
        varify_data(1,i)=lenth*cos(theta);           %the X axis;
        varify_data(2,i)=lenth*sin(theta);           %the Y axis;
        varify_data(3,i)=-2*(-lenth*sin(theta))*rand+2;      %
        varify_label(i)=-1;                          %the expected response; class 1;
    else
        beta=-pi*rand;
        len_2=r+wth*(rand-0.5);
        varify_data(1,i)=len_2*cos(beta)+r;
        varify_data(2,i)=len_2*sin(beta)-dis;
        varify_data(3,i)=(5-len_2*sin(beta))*rand;
        varify_label(i)=1;                   %class 2;
    end
      
end
varify_data = varify_data';
varify_label = varify_label';
pos = find(varify_label == 1);
neg = find(varify_label == -1);
figure;
plot3(varify_data(pos,1),varify_data(pos,2),varify_data(pos,3),'r+');
hold on;
plot3(varify_data(neg,1),varify_data(neg,2),varify_data(neg,3),'b*');
hold off;
title('��֤��');
Kernel = zeros(N,NS);
for j_inside =1:NS
    for i_outside = 1:N
        Kernel(i_outside,j_inside) =exp(-1/(2*sig)*norm(data(i_outside,:)-varify_data(j_inside,:))^2);
    end
end
for k =1:NS
predict(k) = (alphas.*label)'*Kernel(:,k);
end
predict = sign(predict);

% for i = 1:NS
% pre_Li = (alphas.*label)'*(K(data,varify_data,i,sig)) + b;
% pre_Li = sign(pre_Li);
% pre_L(i) = pre_Li;
% end
for i = 1:NS
pre_Li = (alphas.*label)'*(K(data,varify_data,i,sig)) + b;
pre_Li = sign(pre_Li);
pre_L(i) = pre_Li;
    if pre_Li == varify_label(i)
        correct = correct + 1;
    end
end
Accuracy = correct/NS;
fprintf("The accuracy is %.2f%%",Accuracy*100);
pos = find(pre_L == 1);
neg = find(pre_L == -1);
figure;
plot3(varify_data(pos,1),varify_data(pos,2),varify_data(pos,3),'r+');
hold on;
plot3(varify_data(neg,1),varify_data(neg,2),varify_data(neg,3),'b*');
hold off;
title(['Ԥ������׼ȷ��=',num2str(Accuracy*100),'%']);


% �����ֽ��棬�Լ�b��������1�ķֽ���
% hold on
% k = -W(1)/W(2);
% x = -15:0.1:25;
% y = k*x + b;
% plot(x,y,x,y-1,'r--',x,y+1,'r--');
