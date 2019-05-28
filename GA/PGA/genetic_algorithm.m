% Genetic Algorithm 
% C = log2|I + (Es/(M*N0)*(H'*H)|
% population_size: ������Ⱥ��С
% chromosome_size: ����Ⱦɫ�峤��
% generation_size: �����������
% cross_rate: ���뽻�����
% mutate_rate: ����������
% elitism: �����Ƿ�Ӣѡ��
% m: �����Ѹ���
% n: ��������Ӧ��
% p: �����Ѹ�����ֵ�������
% q: �����Ѹ����Ա���ֵ
function [m,n,p] = genetic_algorithm(num_antenna, population_size,...
    chromosome_size, generation_size, cross_rate, mutate_rate,...
    elitism, Es)

global G ;              % ��ǰ��������
global fitness_value;   % ��ǰ����Ӧ�Ⱦ���
global best_fitness;    % �������Ӧֵ
global best_individual; % ��Ѹ���
global best_generation; % ��Ѹ�����ִ�
global good_fitness;    % ���������Ӧ��


fitness_value = zeros(population_size,1);
best_fitness = 0.;
best_generation = 0;

init(population_size, chromosome_size); % ��ʼ��
user = round(3 * rand + 1);                                                     % ��user���û�
H = 1/sqrt(2)*(randn(4,chromosome_size)+1i*randn(4,chromosome_size));           % ��������˹�ֲ�
good_fitness = zeros(1,generation_size);
for G=1:generation_size   
    fitness(num_antenna,user,H,population_size, chromosome_size, Es);           % ������Ӧ�� 
    rank(population_size, chromosome_size);                                     % �Ը��尴��Ӧ�ȴ�С��������
    selection(population_size, chromosome_size, elitism);                       % ѡ�����
    crossover(population_size, chromosome_size, cross_rate);                    % �������
    mutation(population_size, chromosome_size, mutate_rate);                    % �������
end

plotGA(generation_size);% ��ӡ�㷨��������

m = best_individual;    % �����Ѹ���
n = best_fitness;       % ��������Ӧ��
p = best_generation;    % �����Ѹ������ʱ�ĵ�������





