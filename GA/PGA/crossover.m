% 单点交叉操作
% population_size: 种群大小
% chromosome_size: 染色体长度
% cross_rate: 交叉概率

function crossover(population_size, chromosome_size, cross_rate)
global population;
population_new = zeros(population_size, chromosome_size);
mask = zeros(1,chromosome_size);

for i = 1:2:population_size - 5
    for j = 1:chromosome_size
        if rand < cross_rate
            mask(1,j) = 1;
        else
            mask(1,j) = 0;
        end                 
    end                     % 生成掩码
    temp_population = population(i,:) .* mask; 
    population_new(i,:) =  population(i+1,:) .* mask; 
    population_new(i+1,:) = temp_population;
    [~,refill_index] = sort(population(i,:),'descend');
     for l = 1:chromosome_size
        refill_num = 0;
        for m = 1:chromosome_size
           if  population_new(i,m) == chromosome_size - l + 1     %not num 1
               refill_num = 1;
           end
           if refill_num == 1
               break;
           end
        end
        if refill_num == 0
            for k = 1:chromosome_size
                if population_new(i,refill_index(k)) == 0
                   population_new(i,refill_index(k))...
                       = chromosome_size - l + 1;
                   break;
                end
            end
        end
     end
     [~,refill_index] = sort(population(i+1,:),'descend');
     for l = 1:chromosome_size
        refill_num = 0;
        for m = 1:chromosome_size
           if  population_new(i+1,m) == chromosome_size - l + 1     %not num 1
               refill_num = 1;
           end
           if refill_num == 1
               break;
           end
        end
        if refill_num == 0
            for k = 1:chromosome_size
                if population_new(i+1,refill_index(k)) == 0
                   population_new(i+1,refill_index(k))...
                       = chromosome_size - l + 1;
                   break;
                end
            end
        end
     end
end
temp_population = population(population_size,:);
for i = population_size - 3:population_size
    population_new(i,:) = temp_population;
end
population = population_new;

