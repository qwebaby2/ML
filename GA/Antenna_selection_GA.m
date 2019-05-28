clear
clc
%% 生成数据集
%遗传参数
population_size = 120;              % 种群大小
chromosome_size = 32;               % 染色体长度(天线阵列数量)
generation_size = 100;              % 最大迭代次数
cross_rate = 0.6;                   % 交叉概率
mutate_rate = 0.3;                  % 变异概率
elite = zeros(1,chromosome_size);   % 精英
num_antenna = 5;                    % 选择天线数量
sigma = 1;                          % 服从σ的瑞利分布
K = 4;                              % 为4名用户提供服务
P = 1000;                            % 发射总功率
gamma = 300;
user = round(rand * 3 + 1);         % 为第user名用户提供服务
h = raylrnd(sigma,chromosome_size,K);            % 服从瑞利分布的信道矩阵
antenna = zeros(population_size,chromosome_size);% 保存天线选择
%初始化种群
population = zeros(population_size,chromosome_size);
for i=1:population_size
    antenna_index = randperm(32);
    antenna(i,:) = antenna_index;
    for j=1:num_antenna
        % 给population的i行antenna_index(j)列赋值
        population(i,antenna_index(j)) = 1;  
        % 代表第antenna_index(j)个天线被选择
    end
end

%选取SINR为目标函数，信干噪比越大表示该个体适应度越高。
%% 适应度函数
sumh = zeros(chromosome_size,chromosome_size);
fitness_value = zeros(population_size,1);% 所有种群个体适应度初始化为0
for i = 1:K
    sumh = sumh + h(:,i) * h(:,i)';
end

w_ = zeros(chromosome_size,K); %w_是单位波束成形方向矩阵,此时SINR最大
M = zeros(K,K);                %计算功率的中间量
for k = 1:K    
w_(:,k) = (inv(eye(chromosome_size) + (P/(K*sigma))*sumh) * h(:,k)) / ...
    (norm(inv(eye(chromosome_size) + (P/(K*sigma))*sumh) * h(:,k))); %#ok<MINV>
end
for i = 1:K
    for j = 1:K
        if i == j
            M(i,j) = (1/gamma)*(abs(h(:,i)'*w_(:,i)))^2;
        else
            M(i,j) = -(abs(h(:,i)'*w_(:,j)))^2;
        end
    end
end
p = inv(M) * (sigma * ones(K,1)); %#ok<MINV>
w = (sqrt(p) .* (w_'))';

denominator = 0;
%%%%%    SINR = (abs(h(:,user)' * w(:,user))^2) / (denominator + sigma);
for generation = 1:generation_size
    for i = 1:population_size
        denominator = 0;
        for j = 1:K
            if j == user
                continue
            else
                denominator = denominator + (abs(h(:,user)'.* population(i,:) * w(:,j)))^2;
            end
        end
        fitness_value(i) = (abs(h(:,user)'.* population(i,:) * w(:,user))^2)...
            / (denominator + sigma);           
    end
    %% 遍历种群 冒泡排序
    fitness_sum = zeros(population_size,1);
    temp_chromosome = zeros(chromosome_size,1);
    temp_antenna = zeros(1,chromosome_size);
    % 把所有适应度和染色体按适应度从小到大排列
    for i=1:population_size
        min_index = i;
        for j = i+1:population_size
            if fitness_value(j) < fitness_value(min_index)
                min_index = j;
            end
        end
        if min_index ~= i
            % 交换 fitness_value(i) 和 fitness_value(min_index) 的值
            temp = fitness_value(i);
            fitness_value(i) = fitness_value(min_index);
            fitness_value(min_index) = temp;
            % 此时 fitness_value(i) 的适应度在[i,population_size]上最小

            % 交换 population(i) 和 population(min_index) 的染色体串
            for k = 1:chromosome_size
                temp_chromosome(k) = population(i,k);
                population(i,k) = population(min_index,k);
                population(min_index,k) = temp_chromosome(k);
                temp_antenna(k) = antenna(i,k);
                antenna(i,k) = antenna(min_index,k);
                antenna(min_index,k) = temp_antenna(k);
            end
        end
    end
    elite = population(population_size,:);
    % fitness_sum(i) = 前i个个体的适应度之和
    for i=1:population_size
        if i==1
            fitness_sum(i) = fitness_sum(i) + fitness_value(i);    
        else
            fitness_sum(i) = fitness_sum(i-1) + fitness_value(i);
        end
    end
%     %停止迭代条件
%     if generation == 1
%             generation_fitness_avg = fitness_sum(population_size);
%         else
%             lastgeneration_fitness_avg = generation_fitness_avg;
%             generation_fitness_avg = fitness_sum(population_size);
%             dif = abs(generation_fitness_avg - lastgeneration_fitness_avg);
%         if dif < 0.01
%             break;
%         end
%     end
%     
    %% 选择复制
    %锦标赛选择法
    j = 1;
    population_new = zeros(population_size,chromosome_size);
    antenna_new = zeros(population_size,chromosome_size);
    fitness_new = zeros(population_size,1);
    index = randperm(120);
    for i = 1:population_size
        population_new(i,:) = population(index(i),:);
        antenna_new(i,:) = antenna(index(i),:);
        fitness_new(i) = fitness_value(index(i));
    end
    for i = 1:population_size - 1
        if  fitness_new(i) > fitness_new(i+1)
            fitness_value(j) = fitness_new(i);
            antenna(j,:) = antenna_new(i,:);
            population(j,:) = population_new(i,:);
        else
            fitness_value(j) = fitness_new(i+1);
            antenna(j,:) = antenna_new(i+1,:);
            population(j,:) = population_new(i+1,:);
        end
        j = j + 1;
    end
%     for i = 1:population_size
%         if rand < (fitness_value(i) / fitness_value(population_size))
%             population_new(j,:) = population(i,:);
%             antenna_new(j,:) = antenna(i,:);
%             j = j + 1;
%         end
%     end
%     %选择完后随机复制适应度最高样本满足样本数
%     for k = j:population_size
%         %copy_position = round((rand * (j - 1)) + 0.5); 
%         population_new(k,:) = population_new(j - 1,:);
%         antenna_new(k,:) = antenna_new(j - 1,:);
%     end
%     %乱序
%     index = randperm(population_size);
%     for i = 1:population_size
%         population(i,:) =  population_new(index(i),:);
%         antenna(i,:) = antenna_new(index(i),:);
%     end
%     population_new = zeros(population_size,chromosome_size);
%     for i=1:population_size
%         r = rand * fitness_sum(population_size);  
%         % 生成一个随机数，在[0,总适应度]之间
%         first = 1;
%         last = population_size;
%         mid = round((last+first)/2);
%         idx = -1;
% 
%         % 排中法选择个体
%         while (first <= last) && (idx == -1) 
%             if r > fitness_sum(mid)
%                 first = mid;
%             elseif r < fitness_sum(mid)
%                 last = mid;     
%             else
%                 idx = mid;
%                 break;
%             end
%             mid = round((last+first)/2);
%             if (last - first) == 1
%                 idx = last;
%                 break;
%             end
%         end
% 
%        % 产生新一代个体
%        for j=1:chromosome_size
%             population_new(i,j) = population(idx,j);
%        end
%     end
%      population = population_new;

    %% 交叉
    % 步长为2 遍历种群
    index = randperm(120);
    for i = 1:population_size
        population_new(i,:) = population(index(i),:);
        antenna_new(i,:) = antenna(index(i),:);
        fitness_new(i) = fitness_value(index(i));
    end
    population = population_new;
    antenna = antenna_new;
    fitness_value = fitness_new;
    for i=1:2:population_size
        % rand<交叉概率，对两个个体的染色体串进行交叉操作
        if(rand < cross_rate)
            cross_position1 = round((rand * num_antenna) + 0.5);
            cross_position2 = round((rand * num_antenna) + 0.5);
            % 对 两个cross_position处的基因进行交换
            if (population(i,antenna(i,cross_position1)) ==1 && ...
                population(i + 1,antenna(i+1,cross_position1)) == 1) ||...
               (population(i,antenna(i,cross_position2)) ==1 && ...
                population(i + 1,antenna(i+1,cross_position2)) == 1)
                continue
            else
                population(i,antenna(i,cross_position1)) = 0;
                population(i + 1,antenna(i+1,cross_position1)) = 1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                population(i + 1,antenna(i+1,cross_position2)) = 0;
                population(i,antenna(i,cross_position2)) = 1;
                
                temp = antenna(i,cross_position1);
                antenna(i,cross_position1) = antenna(i+1,cross_position2);
                antenna(i+1,cross_position2) = temp;
            end
        end
    end
    %% 变异
    for i=1:population_size
        if rand < mutate_rate
            mutate_position1 = round(rand*(num_antenna + 0.5));
            mutate_position2 = round((chromosome_size - ...
                num_antenna) * rand + 5.5);
            % 变异位置
            if mutate_position1 == 0
                % 若变异位置为0，不变异
                continue;
            end
            population(i,antenna(i,mutate_position1)) = 0;
            population(i,antenna(i,mutate_position2)) = 1;
            temp = antenna(i,mutate_position1);
            antenna(i,mutate_position1) =...
                antenna(i,mutate_position2);
            antenna(i,mutate_position2) = temp;
        end
    end
end
fitness_sum = zeros(population_size,1);
    temp_chromosome = zeros(chromosome_size,1);
    temp_antenna = zeros(1,chromosome_size);
    % 把所有适应度和染色体按适应度从小到大排列
    for i=1:population_size
        min_index = i;
        for j = i+1:population_size
            if fitness_value(j) < fitness_value(min_index)
                min_index = j;
            end
        end
        if min_index ~= i
            % 交换 fitness_value(i) 和 fitness_value(min_index) 的值
            temp = fitness_value(i);
            fitness_value(i) = fitness_value(min_index);
            fitness_value(min_index) = temp;
            % 此时 fitness_value(i) 的适应度在[i,population_size]上最小

            % 交换 population(i) 和 population(min_index) 的染色体串
            for k = 1:chromosome_size
                temp_chromosome(k) = population(i,k);
                population(i,k) = population(min_index,k);
                population(min_index,k) = temp_chromosome(k);
            end
        end
    end
    
denominator1 = 0;
h_index = find(population(population_size,:) == 1);
h_last = zeros(num_antenna,K);
for i = 1:num_antenna
    h_last(i,:) = h(h_index(i),:);
end
w_last = zeros(num_antenna,K);
sumh1 = zeros(num_antenna,num_antenna);
for i = 1:K
    sumh1 = sumh1 + h_last(:,i) * h_last(:,i)';
end
for k = 1:K
w_last(:,k) = (inv(eye(num_antenna) + (P/(K*sigma))*sumh1) * h_last(:,k)) / ...
    (norm(inv(eye(num_antenna) + (P/(K*sigma))*sumh1) * h_last(:,k))); %#ok<MINV>
end
M = zeros(K,K);                %计算功率的中间量
for i = 1:K
    for j = 1:K
        if i == j
            M(i,j) = (1/gamma)*(abs(h_last(:,i)'*w_last(:,i)))^2;
        else
            M(i,j) = -(abs(h_last(:,i)'*w_last(:,j)))^2;
        end
    end
end
p = inv(M) * (sigma * ones(K,1)); %#ok<MINV>
w_last = (sqrt(p) .* (w_last'))';
for j = 1:K
    if j == user
        continue
    else
       denominator1 = denominator1 + abs(h_last(:,user)' * w_last(:,j))^2;
    end
end
fitness_value_last = (abs(h_last(:,user)' * w_last(:,user))^2)...
                     / (denominator1 + 1);           
fprintf('选择的天线下标为：');
h_index %#ok<NOPTS>
fprintf('此时信干噪比为：(dB)');
10 * log10(fitness_value_last) %#ok<NOPTS>
