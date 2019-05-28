function [J,Ej] = select(i,data,num_data,alphas,label,b,C,Ei,choose,sig)
maxDeltaE = 0;maxJ = -1;
if choose == 1 %ȫ����---���ѡ��alphas
    temp = 1;
        while temp
            j = randi(num_data,1);
            temp = (j == i);
        end
    J = j;
    Ej = kernel(data,alphas,label,b,J,sig);
else %���ֱ���--����ʽ��ѡ��alphas
    index = find(alphas>0 & alphas < C);
    for k = 1:length(index)
        if i == index(k)
            continue;
        end
        temp_e = kernel(data,alphas,label,b,k,sig);
        deltaE = abs(Ei - temp_e); %ѡ����Ei�������alphas
        if deltaE > maxDeltaE
            maxJ = k;
            maxDeltaE = deltaE;
            Ej = temp_e;
        end
    end
    J = maxJ;
end