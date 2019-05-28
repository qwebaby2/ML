% Genetic Algorithm 
% C = log2|I + (Es/(M*N0)*(H'*H)|
% population_size: 输入种群大小
% chromosome_size: 输入染色体长度
% generation_size: 输入迭代次数
% cross_rate: 输入交叉概率
% mutate_rate: 输入变异概率
% elitism: 输入是否精英选择
% m: 输出最佳个体
% n: 输出最佳适应度
% p: 输出最佳个体出现迭代次数
% q: 输出最佳个体自变量值
function [m,n,p] = genetic_algorithm(num_antenna, population_size,...
    chromosome_size, generation_size, cross_rate, mutate_rate,...
    elitism, Es)

global G ;              % 当前迭代次数
global fitness_value;   % 当前代适应度矩阵
global best_fitness;    % 总最佳适应值
global best_individual; % 最佳个体
global best_generation; % 最佳个体出现代
global good_fitness;    % 历代最佳适应度


fitness_value = zeros(population_size,1);
best_fitness = 0.;
best_generation = 0;

init(population_size, chromosome_size); % 初始化
user = round(3 * rand + 1);                                                     % 第user名用户
H = 1/sqrt(2)*(randn(4,chromosome_size)+1i*randn(4,chromosome_size));           % 产生复高斯分布
good_fitness = zeros(1,generation_size);
for G=1:generation_size   
    fitness(num_antenna,user,H,population_size, chromosome_size, Es);           % 计算适应度 
    rank(population_size, chromosome_size);                                     % 对个体按适应度大小进行排序
    selection(population_size, chromosome_size, elitism);                       % 选择操作
    crossover(population_size, chromosome_size, cross_rate);                    % 交叉操作
    mutation(population_size, chromosome_size, mutate_rate);                    % 变异操作
end

plotGA(generation_size);% 打印算法迭代过程

m = best_individual;    % 获得最佳个体
n = best_fitness;       % 获得最佳适应度
p = best_generation;    % 获得最佳个体出现时的迭代次数





