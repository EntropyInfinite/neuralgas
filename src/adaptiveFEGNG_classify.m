function classification_rate  = adaptiveFEGNG_classify(nodes, edges, node_classes, node_lambdas, point_coverages, testData, testLabels)
%AFEGNG Multivariate classification rate
%   Computes the multivariate classification rate for neural gas on the
%   testing dataset provided

dataLength = size(testData, 2);
producedLabels = zeros(size(testLabels));
bida = 0;

for n=1:dataLength
    %locate voter nodes (which have this point inside their lambda-radius)
    distances = nodes-repmat(testData(:,n),1,size(nodes,2));
    norms = sqrt(sum(distances.^2,1));
    voter_nodes_indices = find(norms<=node_lambdas);
    
    %handle the case when there are no voters (outlier point)
    if isempty(voter_nodes_indices)
        [nsorted, isorted] = sort(norms);
        ClassIdents = unique(node_classes);
        NumOfClasses = numel(ClassIdents);
        classprobs = zeros(1,NumOfClasses);
        for k=1:NumOfClasses
            ind = find(node_classes(isorted)==ClassIdents(k),1);
            lmb = node_lambdas(isorted);
            classprobs(k) = nsorted(ind)-lmb(ind);
        end
        [value, index] = min(classprobs);
        producedLabels(n) = ClassIdents(index);
    else
        distance_ratings = norms(voter_nodes_indices)/sum(norms(voter_nodes_indices));
        lambda_ratings = node_lambdas(voter_nodes_indices)/sum(node_lambdas(voter_nodes_indices));
        count_ratings = point_coverages(voter_nodes_indices)/sum(point_coverages(voter_nodes_indices));
        [value, index] = min(distance_ratings+lambda_ratings-count_ratings);
        producedLabels(n) = node_classes(voter_nodes_indices(index));
    end
end

comparison = producedLabels == testLabels;
classification_rate = double(sum(comparison))/dataLength;

end

