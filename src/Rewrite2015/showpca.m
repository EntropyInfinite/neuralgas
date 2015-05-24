function [] = showpca( input, inputL )

[coeff,score,latent,tsquared,explained,mu] = pca(input);
display(explained)
gscatter(score(:,1),score(:,2),inputL)
titlestr = ['Total variance explained ' num2str(explained(1)+explained(2)) '%' ];
title(titlestr)
end

