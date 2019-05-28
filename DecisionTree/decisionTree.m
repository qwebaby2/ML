function decisionTreeModel=decisionTree(data,label,propertyName,delta)

global Node;

Node = struct('fatherNodeName',[],'EdgeProperty',[],'NodeName',[]);
BuildTree('root','Stem',data,label,propertyName,delta);
Node(1) = [];
model.Node = Node;
decisionTreeModel = model;
