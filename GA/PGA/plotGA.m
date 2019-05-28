%打印算法迭代过程
%generation_size: 迭代次数

function plotGA(generation_size)
global good_fitness;
x = 1:1:generation_size;
y = good_fitness;
plot(x,y)
