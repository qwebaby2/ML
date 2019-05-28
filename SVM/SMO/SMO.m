%% 开始前的准备
clear
clc
%% 读取数据并整理
load('Optimized_SMO_data.mat');
y(y == 0) = -1;
data = [X,y]'; % 数据集
train_data = data(1:end-1,:)';
label = data(end,:)'; % 标签
[num_data,d] = size(train_data); % 样本总数
data = train_data; % 特征向量集合
%% 定义向量机参数
alphas = zeros(num_data,1); % 拉格朗日乘数个数等于样本个数
% 偏置项
b = 0;
% 松弛变量的最大值 i.e. ksi <= C
C = 0.6;
% 迭代计数变量和最大迭代次数
iter = 0;
max_iter = 40;
% 迭代过程
while iter < max_iter % 设置循环条件
    alpha_change = 0; % 每次迭代初始化乘子改变情况
    for i = 1:num_data % 第一个样本i的迭代
        %样本i的输出目标值
        predicted_label_i = (alphas.*label)'*(data*data(i,:)') + b;
        %样本i误差
        E_i = predicted_label_i - label(i);
        % 满足KKT条件
        if (label(i)*E_i < -0.001 && alphas(i) < C) || (label(i)*E_i > 0.001 && alphas(i) > 0) %epsilon为精度
           % 选择一个和 i 不相同的待改变的alpha(j)
            j = randi(num_data,1);  
            if j == i % 将j变为和i不同的值
                temp = 1;
                while temp % 第一次肯定会执行，直到j和i不同时跳出循环
                    j = randi(num_data,1);
                    if j ~= i
                        temp = 0;
                    end
                end
            end
            % 样本j的目标输出值
            predicted_label_j = (alphas.*label)'*(data*data(j,:)') + b;
            %样本j误差
            E_j = predicted_label_j - label(j);
            %更新上下限
            if label(i) ~= label(j) %类标签不同
                L = max(0,alphas(j) - alphas(i));
                H = min(C,C + alphas(j) - alphas(i));
            else % 类标签相同
                L = max(0,alphas(j) + alphas(i) - C);
                H = min(C,alphas(j) + alphas(i));
            end
            if L == H  %上下限一样结束本次j值的执行，选新的j值进行下一次迭代
                continue;
            end
            %计算eta
            eta = -2*data(i,:)*data(j,:)' + data(i,:)*data(i,:)' + ...
                  data(j,:)*data(j,:)';
            %保存旧值
            alphas_i_old = alphas(i);
            alphas_j_old = alphas(j);
            %更新alpha(j)
            alphas(j) = alphas(j) + label(j)*(E_i-E_j)/eta;
            %限制范围
            if alphas(j) > H
                alphas(j) = H;
            elseif alphas(j) < L
                    alphas(j) = L;
            end
            %如果alpha(j)没怎么改变，结束本次j的执行，选新的j值进行下一次迭代
            if abs(alphas(j) - alphas_j_old) < 1e-4
                continue;
            end
            %更新alpha(i)
            alphas(i) = alphas(i) + label(i)*label(j)*(alphas_j_old-alphas(j));
            %更新系数b
            b1 = b - E_i - label(i)*(alphas(i) - alphas_i_old)*data(i,:)*data(i,:)' - ...
                label(j)*(alphas(j) - alphas_j_old)*data(i,:)*data(j,:)';
            b2 = b - E_j - label(i)*(alphas(i) - alphas_i_old)*data(i,:)*data(j,:)' - ...
                label(j)*(alphas(j) - alphas_j_old)*data(j,:)*data(j,:)';
            %b的几种选择机制
            if alphas(i) > 0 && alphas(i) < C
                b = b1;
            elseif alphas(j) > 0 && alphas(j) < C
                b = b2;
            else
                b = (b1 + b2)/2;
            end
            %确定更新了，记录一次
            alpha_change = alpha_change + 1;
        end
    end
    % 没有实行alpha更新，迭代加1
    if alpha_change == 0
        iter = iter + 1;
    %实行了交换，迭代清0
    else
        iter = 0;
    end
    disp(['iter: ',num2str(iter)]);
end
%% 计算权值w
w = (alphas.*label)'*data;
%记录支持向量位置
index_sv = find(alphas ~= 0);
%计算预测结果
predict = (alphas.*label)'*(data*data') + b;
predict = sign(predict);
%% 显示结果
figure;
neg = find(predict == -1);
data_neg = (data(neg,:))';
plot(data_neg(1,:),data_neg(2,:),'*r'); % 负类
hold on
pos = find(predict == 1);
data_pos = (data(pos,:))';
plot(data_pos(1,:),data_pos(2,:),'+g'); % 正类
hold on
data_sv = (data(index_sv,:))'; % 支持向量
plot(data_sv(1,:),data_sv(2,:),'p','LineWidth',2);
% 画出分界面，以及b上下正负1的分界面
hold on
k = -w(1)/w(2);
b = -b/w(2);
x = 0:0.1:5;
y = k*x + b;
plot(x,y,'k',x,y-1,'b--',x,y+1,'b--');
title(['C = ',num2str(C)]);