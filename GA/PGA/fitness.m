% ������Ⱥ������Ӧ�ȣ��Բ�ͬ���Ż�Ŀ�꣬�޸�����ĺ���
% population_size: ��Ⱥ��С
% chromosome_size: Ⱦɫ�峤��
% ѡ��C = log2|I + (Es/(M*N0)*(H'*H)| ��λ������������޴�������������
function fitness(num_antenna, user,H,population_size, chromosome_size, Es)
global fitness_value;
global population;
global H_;

% ������Ⱥ������Ӧ�ȳ�ʼ��Ϊ0
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

