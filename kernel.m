function Ek = kernel(data1,alphas,label,b,k,sig)
pre_Li = (alphas.*label)'*(K(data1,data1,k,sig)) + b;
%pre_Li = (alphas.*label)'*(data*data(k,:)') + b;
Ek = pre_Li - label(k);