% ��ʼ����Ⱥ
% population_size: ��Ⱥ��С
% chromosome_size: Ⱦɫ�峤��

function init(population_size, chromosome_size)
global population;
for i=1:population_size
        population(i,:) = randperm(chromosome_size);  % rand����(0,1)֮����������round()���������뺯��
end
clear i;



