function classification_rate  = adaptiveFEGNG_binclass(nodes, edges, node_classes, node_lambdas, point_coverages, testData, testLabels)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dataLength = size(testData, 2);

winner_indices = knnsearch(nodes', testData');
comparison = node_classes(winner_indices)' == testLabels;
classification_rate = double(sum(comparison))/dataLength;

end

