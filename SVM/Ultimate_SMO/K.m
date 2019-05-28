function f=K(data1,data2,z,sig)
%sigma=0.25;
num_data1 = length(data1);
f = zeros(num_data1,1);
for i = 1:num_data1
    f(i) = exp(-1/(2*sig)*norm(data1(i,:)-data2(z,:))^2);
end
