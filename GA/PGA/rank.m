% 对个体按适应度大小进行排序，并且保存最佳个体
% population_size: 种群大小
% chromosome_size: 染色体长度

function rank(population_size, chromosome_size)
global fitness_value;   % 种群适应度
global best_fitness;
global best_individual;
global best_generation;
global population;
global good_fitness;
global G;


temp_chromosome(chromosome_size)=0;

% 遍历种群 
% 冒泡排序
% 最后population(i)的适应度递增
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
good_fitness(1,G) = fitness_value(population_size);
% 更新最大适应度和对应的迭代次数，保存最佳个体(最佳个体的适应度最大)
if fitness_value(population_size) > best_fitness
    best_fitness = fitness_value(population_size);
    best_generation = G;
    for j=1:chromosome_size
        best_individual(j) = population(population_size,j);
    end
end

