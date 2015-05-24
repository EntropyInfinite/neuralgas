function cv_result = crossval_AGNHc(Set, Labels, folds)
    close all
    
    inds = randperm(size(Set,1));
    DSample = Set(inds,:);
    CL = Labels(inds);
    
    foldsize = floor(size(Set,1)/folds);
    foldind = [0 (1:folds-1)*foldsize size(Set,1)];
    
    for nfold=1:folds
        ind1 = foldind(nfold);
        ind2 = foldind(nfold+1);
        testSet = DSample(ind1+1:ind2,:);
        testLabels = CL(ind1+1:ind2,:);
        trainSet = [DSample(1:ind1,:); DSample(ind2+1:end,:)];
        trainLabels = [CL(1:ind1,:); CL(ind2+1:end,:)];
        cv_result.fold_results(nfold) = AGNHc_classtest(trainSet, trainLabels, testSet, testLabels);
        close all
    end
    
    cv_result.fitdist = fitdist([cv_result.fold_results.('accuracy')]', 'normal');
    cv_result.conf_ints = paramci(cv_result.fitdist);
end

