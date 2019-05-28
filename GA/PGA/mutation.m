% 单点变异操作
% population_size: 种群大小
% chromosome_size: 染色体长度
% mutate_rate: 变异概率

function mutation(population_size, chromosome_size, mutate_rate)
global population;
mask = zeros(1,chromosome_size);
for i=1:population_size - 1
    for j = 1:chromosome_size
        if rand < mutate_rate
            mask(j) = 1;
        end
    end
    for j = 1:chromosome_size
        if mask(j) == 1
            comp = 1;
            while comp
                index = round((31 * rand) + 1);
                if index ~= j
                    comp = 0;
                end
            end
            temp = population(i,j);
            population(i,j) = population(i,index);
            population(i,index) = temp;
        end
    end
end



