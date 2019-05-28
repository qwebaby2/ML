% ������ѡ��
% population_size: ��Ⱥ��С
% chromosome_size: Ⱦɫ�峤��
% elitism: �Ƿ�Ӣѡ��

function selection(population_size, chromosome_size, elitism)

 global population;          % ǰ����Ⱥ
 global fitness_value;       % ǰ����Ӧ��
% global population_new;      % ��һ����Ⱥ
% global fitness_value_new;   % ��һ����Ӧ��
population_new = zeros(population_size, chromosome_size);
fitness_value_new = zeros(population_size,1);
if elitism == 0
    for i=1:population_size
       s = round((119 * rand(1,3)) + 1);    % ���ѡ��3��Ⱦɫ�����ѡ���
       s = [fitness_value(s(1)),fitness_value(s(2)),fitness_value(s(3))];
       chosen = find(fitness_value == max(s));
       population_new(i,:) = population(chosen,:); 
       fitness_value_new(i,:) = fitness_value(chosen,:);
    end
else
    for i=1:population_size - 5
       s = round((119 * rand(1,3)) + 1);    % ���ѡ��3��Ⱦɫ�����ѡ���
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


 

