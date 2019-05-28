function label = decisionTreeTest(decisionTreeModel,sampleSet,propertyName)

lengthSample = size(sampleSet,1);
label = -1 * ones(lengthSample,1);
Nodes = decisionTreeModel.Node;
Nodes(2).level = 1;
Nodes(3).level = 2;
Nodes(4).level = 2;
Nodes(5).level = 1;
Nodes(6).level = 2;
Nodes(7).level = 2;
Nodes(8).level = 2;
for sampleIndex = 1:lengthSample
    sample = sampleSet(sampleIndex,:);
    rootNode = Nodes(1);
    head = rootNode.NodeName;
    index = GetFeatureNum(propertyName,head);
    edge = sample(index);
    k = 1;
    level = 1;
    while k < length(Nodes)
        k = k + 1;
        if Nodes(k).level == level
            if strcmp(Nodes(k).fatherNodeName,head)
                if Nodes(k).EdgeProperty == edge
                    if Nodes(k).NodeName < 10
                        label(sampleIndex) = Nodes(k).NodeName;
                        break;
                    else
                        head=Nodes(k).NodeName;
                        index=GetFeatureNum(propertyName,head);
                        edge=sample(index);
                        level=level+1;
                    end
                end
            end
        end
    end
end
