% �Ը��尴��Ӧ�ȴ�С�������򣬲��ұ�����Ѹ���
% population_size: ��Ⱥ��С
% chromosome_size: Ⱦɫ�峤��

function rank(population_size, chromosome_size)
global fitness_value;   % ��Ⱥ��Ӧ��
global best_fitness;
global best_individual;
global best_generation;
global population;
global good_fitness;
global G;


temp_chromosome(chromosome_size)=0;

% ������Ⱥ 
% ð������
% ���population(i)����Ӧ�ȵ���
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
good_fitness(1,G) = fitness_value(population_size);
% ���������Ӧ�ȺͶ�Ӧ�ĵ���������������Ѹ���(��Ѹ������Ӧ�����)
if fitness_value(population_size) > best_fitness
    best_fitness = fitness_value(population_size);
    best_generation = G;
    for j=1:chromosome_size
        best_individual(j) = population(population_size,j);
    end
end

