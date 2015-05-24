function scores = AGNHc_classtest(trainSet, trainLabels, testSet, testLabels)
    close all
    lb = unique(trainLabels);
    figure
    a = AGNHcE(trainSet,trainLabels,size(trainSet,1));
    scores.main = a;
    predicted = AGNHc_classifyL(a,testSet);
    subplot(2,1,2)
    bar(testLabels)
    [~,foo] = max(predicted,[],2);
    cleanpred = lb(foo);
    scores.accuracy = sum(cleanpred==testLabels)/length(testLabels);
    [scores.confmatrix, scores.confmatOrder] = confusionmat(testLabels, cleanpred);
    if length(lb)>2
        return
    end
    %[X,Y,~,AUC] = perfcurve(testLabels,predicted(:,1),lb(1));
    %figure
    %plot(X,Y)
    %scores.AUC = AUC;
end

