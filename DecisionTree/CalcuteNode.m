function [NodeIndex,BuildNode]=CalcuteNode(data,label,delta)

LargeEntropy = CEntropy(label);
[m,n] = size(data);
EntropyGain = LargeEntropy*ones(1,n);
BuildNode = true;
for i = 1:n
    pData = data(:,i);
    itemList = unique(pData);
    for j = 1:length(itemList)
        itemIndex = find(pData == itemList(j));
        EntropyGain(i) = EntropyGain(i) - (length(itemIndex)) / m * ...
            CEntropy(label(itemIndex));
    end
    % 此处运行则为增益率，注释掉则为增益
    % EntropyGain(i) = EntropyGain(i)/CEntropy(pData); 
end
[maxGainEntropy,NodeIndex] = max(EntropyGain);
if maxGainEntropy < delta
    BuildNode = false;
end
