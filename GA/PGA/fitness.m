% 计算种群个体适应度，对不同的优化目标，修改下面的函数
% population_size: 种群大小
% chromosome_size: 染色体长度
% 选择C = log2|I + (Es/(M*N0)*(H'*H)| 单位带宽内允许的无错误的最大传输码率
function fitness(num_antenna, user,H,population_size, chromosome_size, Es)
global fitness_value;
global population;
global H_;

% 所有种群个体适应度初始化为0
fitness_value = zeros(population_size,1);    
index = zeros(population_size,chromosome_size);
H_ = zeros(1,num_antenna);
for i = 1:population_size
    [~,index(i,:)] = sort(population(i,:),'descend');
end


% C = log2|I + (Es / (M*N0)*(H' * H)|;
for i = 1:population_size
    for j = 1:num_antenna
        H_(j) = H(user,index(i,j));
    end
    fitness_value(i) = real(log2(det(eye(num_antenna) + (Es/(num_antenna...
        *0.75)*(H_' * H_)))));
end

