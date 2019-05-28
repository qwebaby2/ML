clear
close
clc
%% 
data = load('C:\Users\dell\Desktop\code\input.csv');
[row, col] = size(data);
features = data(:,1:col-1);
label = data(:,col);
for ii = 1:row
    if label(ii) == 10
        label(ii) = -1;
    else
        label(ii) = 1;
    end
end
[features_norm,PS] = mapminmax(features',-1,1);
features_norm = features_norm';
[X_train, y_train, X_test, y_test] = split_train(features_norm, label, 2, 0.7);
model = svmtrain(y_train, X_train, '-c 0.5 -g 1 -v 5 -h 0');
[pred, acc, ~] = svmpredict(y_test, X_test, model);
[bestacc,bestc,bestg] = SVMcg(y_train,X_train,-5,8,-8,8,10,0.5,0.5,0.9);
cmd = ['-c',num2str(bestc),'-g',num2str(bestg),'-w1 2 -w0 0.5'];
