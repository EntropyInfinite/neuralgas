function [sdistances, indices] = findNodeDistances(nodes) %#codegen

NumOfNodes = size(nodes,2);

distances = zeros(NumOfNodes);
for i=1:NumOfNodes
    for j=1:NumOfNodes
       distances(i,j) = norm(nodes(:,j) - nodes(:,i)); 
    end
end

[sdistances, indices] = sort(distances,2);

%sdistances = sdistances(2:end, 2:end);
%indices = indices(2:end, 2:end);

% MEX Code Generation:

% mexcfg = coder.config('mex');
% mexcfg.DynamicMemoryAllocation = 'AllVariableSizeArrays'; 
% codegen -config mexcfg findTwoNearest -args {coder.typeof(In(:,n),[Inf 1]),coder.typeof(nodes,[Inf Inf])}