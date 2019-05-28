
function main()
clear
clc
elitism = 1;                % 选择精英操作
population_size = 120;      % 种群大小
chromosome_size = 32;       % 染色体长度
generation_size = 150;       % 最大迭代次数
cross_rate = 0.6;           % 交叉概率
mutate_rate = 0.09;         % 变异概率
num_antenna = 5;            % 选择天线数目
Es = 100;                   % 传输功率
[best_individual,best_fitness,iterations] = genetic_algorithm...
    (num_antenna, population_size, chromosome_size, generation_size,...
    cross_rate, mutate_rate,elitism, Es);
[~,index] = sort(best_individual,'descend');

disp 最优天线选择:
sort(index(1:5))
disp 最优适应度:
best_fitness
disp 达到最优结果的迭代次数:
iterations

clear