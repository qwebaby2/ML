
function main()
clear
clc
elitism = 1;                % ѡ��Ӣ����
population_size = 120;      % ��Ⱥ��С
chromosome_size = 32;       % Ⱦɫ�峤��
generation_size = 150;       % ����������
cross_rate = 0.6;           % �������
mutate_rate = 0.09;         % �������
num_antenna = 5;            % ѡ��������Ŀ
Es = 100;                   % ���书��
[best_individual,best_fitness,iterations] = genetic_algorithm...
    (num_antenna, population_size, chromosome_size, generation_size,...
    cross_rate, mutate_rate,elitism, Es);
[~,index] = sort(best_individual,'descend');

disp ��������ѡ��:
sort(index(1:5))
disp ������Ӧ��:
best_fitness
disp �ﵽ���Ž���ĵ�������:
iterations

clear