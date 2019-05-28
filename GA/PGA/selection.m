% 锦标赛选择法
% population_size: 种群大小
% chromosome_size: 染色体长度
% elitism: 是否精英选择

function selection(population_size, chromosome_size, elitism)

 global population;          % 前代种群
 global fitness_value;       % 前代适应度
% global population_new;      % 新一代种群
% global fitness_value_new;   % 新一代适应度
population_new = zeros(population_size, chromosome_size);
fitness_value_new = zeros(population_size,1);
if elitism == 0
    for i=1:population_size
       s = round((119 * rand(1,3)) + 1);    % 随机选择3个染色体进入选择池
       s = [fitness_value(s(1)),fitness_value(s(2)),fitness_value(s(3))];
       chosen = find(fitness_value == max(s));
       population_new(i,:) = population(chosen,:); 
       fitness_value_new(i,:) = fitness_value(chosen,:);
    end
else
    for i=1:population_size - 5
       s = round((119 * rand(1,3)) + 1);    % 随机选择3个染色体进入选择池
       s = [fitness_value(s(1)),fitness_value(s(2)),fitness_value(s(3))];
       chosen = find(fitness_value == max(s));
       [h1,~] = size(chosen);
       if h1 ~= 1
           chosen = chosen(h1,1);           
       end
       population_new(i,:) = population(chosen,:); 
       fitness_value_new(i,:) = fitness_value(chosen,:);
    end
    for i=population_size - 4:population_size 
        population_new(i,:) = population(population_size,:); 
        fitness_value_new(i,:) = fitness_value(population_size,:);
    end
%     diso = randperm(population_size);
%     popualtion_diso = zeros(population_size,chromosome_size);
%     fitness_value_diso = zeros(population_size,1);
%     for i = 1:population_size
%        popualtion_diso(i,:) = population_new(diso(i),:);
%        fitness_value_diso(i,:) = fitness_value_new(diso(i),:);
%     end
%     population_new = popualtion_diso;
%     fitness_value_new = fitness_value_diso;
end
population = population_new;
fitness_value = fitness_value_new; 


 

