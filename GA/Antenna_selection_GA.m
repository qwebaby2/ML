clear
clc
%% �������ݼ�
%�Ŵ�����
population_size = 120;              % ��Ⱥ��С
chromosome_size = 32;               % Ⱦɫ�峤��(������������)
generation_size = 100;              % ����������
cross_rate = 0.6;                   % �������
mutate_rate = 0.3;                  % �������
elite = zeros(1,chromosome_size);   % ��Ӣ
num_antenna = 5;                    % ѡ����������
sigma = 1;                          % ���Ӧҵ������ֲ�
K = 4;                              % Ϊ4���û��ṩ����
P = 1000;                            % �����ܹ���
gamma = 300;
user = round(rand * 3 + 1);         % Ϊ��user���û��ṩ����
h = raylrnd(sigma,chromosome_size,K);            % ���������ֲ����ŵ�����
antenna = zeros(population_size,chromosome_size);% ��������ѡ��
%��ʼ����Ⱥ
population = zeros(population_size,chromosome_size);
for i=1:population_size
    antenna_index = randperm(32);
    antenna(i,:) = antenna_index;
    for j=1:num_antenna
        % ��population��i��antenna_index(j)�и�ֵ
        population(i,antenna_index(j)) = 1;  
        % �����antenna_index(j)�����߱�ѡ��
    end
end

%ѡȡSINRΪĿ�꺯�����Ÿ����Խ���ʾ�ø�����Ӧ��Խ�ߡ�
%% ��Ӧ�Ⱥ���
sumh = zeros(chromosome_size,chromosome_size);
fitness_value = zeros(population_size,1);% ������Ⱥ������Ӧ�ȳ�ʼ��Ϊ0
for i = 1:K
    sumh = sumh + h(:,i) * h(:,i)';
end

w_ = zeros(chromosome_size,K); %w_�ǵ�λ�������η������,��ʱSINR���
M = zeros(K,K);                %���㹦�ʵ��м���
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
    %% ������Ⱥ ð������
    fitness_sum = zeros(population_size,1);
    temp_chromosome = zeros(chromosome_size,1);
    temp_antenna = zeros(1,chromosome_size);
    % ��������Ӧ�Ⱥ�Ⱦɫ�尴��Ӧ�ȴ�С��������
    for i=1:population_size
        min_index = i;
        for j = i+1:population_size
            if fitness_value(j) < fitness_value(min_index)
                min_index = j;
            end
        end
        if min_index ~= i
            % ���� fitness_value(i) �� fitness_value(min_index) ��ֵ
            temp = fitness_value(i);
            fitness_value(i) = fitness_value(min_index);
            fitness_value(min_index) = temp;
            % ��ʱ fitness_value(i) ����Ӧ����[i,population_size]����С

            % ���� population(i) �� population(min_index) ��Ⱦɫ�崮
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
    % fitness_sum(i) = ǰi���������Ӧ��֮��
    for i=1:population_size
        if i==1
            fitness_sum(i) = fitness_sum(i) + fitness_value(i);    
        else
            fitness_sum(i) = fitness_sum(i-1) + fitness_value(i);
        end
    end
%     %ֹͣ��������
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
    %% ѡ����
    %������ѡ��
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
%     %ѡ��������������Ӧ�������������������
%     for k = j:population_size
%         %copy_position = round((rand * (j - 1)) + 0.5); 
%         population_new(k,:) = population_new(j - 1,:);
%         antenna_new(k,:) = antenna_new(j - 1,:);
%     end
%     %����
%     index = randperm(population_size);
%     for i = 1:population_size
%         population(i,:) =  population_new(index(i),:);
%         antenna(i,:) = antenna_new(index(i),:);
%     end
%     population_new = zeros(population_size,chromosome_size);
%     for i=1:population_size
%         r = rand * fitness_sum(population_size);  
%         % ����һ�����������[0,����Ӧ��]֮��
%         first = 1;
%         last = population_size;
%         mid = round((last+first)/2);
%         idx = -1;
% 
%         % ���з�ѡ�����
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
%        % ������һ������
%        for j=1:chromosome_size
%             population_new(i,j) = population(idx,j);
%        end
%     end
%      population = population_new;

    %% ����
    % ����Ϊ2 ������Ⱥ
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
        % rand<������ʣ������������Ⱦɫ�崮���н������
        if(rand < cross_rate)
            cross_position1 = round((rand * num_antenna) + 0.5);
            cross_position2 = round((rand * num_antenna) + 0.5);
            % �� ����cross_position���Ļ�����н���
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
    %% ����
    for i=1:population_size
        if rand < mutate_rate
            mutate_position1 = round(rand*(num_antenna + 0.5));
            mutate_position2 = round((chromosome_size - ...
                num_antenna) * rand + 5.5);
            % ����λ��
            if mutate_position1 == 0
                % ������λ��Ϊ0��������
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
    % ��������Ӧ�Ⱥ�Ⱦɫ�尴��Ӧ�ȴ�С��������
    for i=1:population_size
        min_index = i;
        for j = i+1:population_size
            if fitness_value(j) < fitness_value(min_index)
                min_index = j;
            end
        end
        if min_index ~= i
            % ���� fitness_value(i) �� fitness_value(min_index) ��ֵ
            temp = fitness_value(i);
            fitness_value(i) = fitness_value(min_index);
            fitness_value(min_index) = temp;
            % ��ʱ fitness_value(i) ����Ӧ����[i,population_size]����С

            % ���� population(i) �� population(min_index) ��Ⱦɫ�崮
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
M = zeros(K,K);                %���㹦�ʵ��м���
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
fprintf('ѡ��������±�Ϊ��');
h_index %#ok<NOPTS>
fprintf('��ʱ�Ÿ����Ϊ��(dB)');
10 * log10(fitness_value_last) %#ok<NOPTS>
