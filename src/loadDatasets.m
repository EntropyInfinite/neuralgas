liverTrain = liver(1:250,1:6)';
liverTrainL = liver(1:250,7);
liverTest = liver(251:345,1:6)';
liverTestL = liver(251:345,7);
failscript
[ nodes, edges, node_classes ] = construct_fastGNG_classifier(liverTrain, liverTrainL, 1, 250, 1, 10.2);
[ nodes, edges, node_classes ] = construct_fastGNG_classifier(cancerTrain, cancerTrainL, 1, 500, 1, 4.7);
fastGNG_classify(nodes, edges, node_classes, cancerTest, cancerTestL);
clear
ionTrain = ion(1:260,1:34)';
ionTrainL = ion(1:260,35);
ionTest = ion(261:351,1:34)';
ionTestL = ion(261:351,35);
failscript
[ nodes, edges, node_classes ] = construct_fastGNG_classifier(ionTrain, ionTrainL, 1, 260, 1, 1.09);
ion = ion(randperm(size(ion,1)),:);
ionTrain = ion(1:260,1:34)';
ionTrainL = ion(1:260,35);
ionTestL = ion(261:351,35);
ionTest = ion(261:351,1:34)';
failscript
liver = liver(randperm(size(liver,1)),:);
liverTrain = liver(1:250,1:6)';
liverTrainL = liver(1:250,7);
liverTest = liver(251:345,1:6)';
liverTestL = liver(251:345,7);
failscript
heart = heart(randperm(size(heart,1)),:);
heartTrain = heart(1:220,1:13)';
heartTest = heart(221:297,1:13)';
heartTrainL = heart(1:220,14)>0;
heartTestL = heart(221:297,14)>0;
heartTestL = +heartTestL;
heartTrainL = +heartTrainL;

hepa = hepa(randperm(size(hepa,1)),:);
hepaTrain =hepa(1:110,2:20)';
hepaTest =hepa(111:155,2:20)';
hepaTrainL =hepa(1:110,1);
hepaTestL =hepa(111:155,1);


liver = liver(randperm(size(liver,1)),:);
liverTrain = liver(1:250,1:6)';
liverTrainL = liver(1:250,7);
liverTest = liver(251:345,1:6)';
liverTestL = liver(251:345,7);
failscript
[ nodes, edges, node_classes ] = construct_EfastGNG_classifier(balanceTrain, balanceTrainL, 1, 460, lambda);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier(ionTrain, ionTrainL, 1, 260, 1);
fastGNG_classify(nodes, edges, node_classes, ionTest, ionTestL);
failscript
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier(ionTrain, ionTrainL, 1, 260, 0);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier(ionTrain, ionTrainL, 1, 260, 0.9);
construct_fastEGNG_classifier
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier(ionTrain, ionTrainL, 1, 260, 0.9);
fastGNG_classify(nodes, edges, node_classes, ionTest, ionTestL);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier(ionTrain, ionTrainL, 1, 260, 0.5);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier(ionTrain, ionTrainL, 1, 260, 0.8);
fastGNG_classify(nodes, edges, node_classes, ionTest, ionTestL);
failscript
cancer = cancer(randperm(size(cancer,1)),:);
cancerTrain = cancer(1:500,1:9)';
cancerTrainL = cancer(1:500,10);
cancerTest = cancer(501:683,1:9)';
cancerTestL = cancer(501:683,10);
failscript
balance = balance(randperm(size(balance,1)),:);
failscript
balanceTrain = balancescale(1:460,2:5)';
balanceTrainLabel = balancescale(1:460,1);
balanceTest = balancescale(461:625,2:5)';
balanceTestLabel = balancescale(461:625,1);
failscript
newthyroid = newthyroid(randperm(size(newthyroid,1)),:);
newthyroid1 = mapstd(newthyroid(:,2:6)')';
thyrTrain = newthyroid1(1:160,:)';
thyrTest = newthyroid1(161:215,:)';
thyrTrainL = newthyroid(1:160,1);
thyrTestL = newthyroid(161:215,1);
failscript
diab = diab(randperm(size(diab,1)),:);
diab1 = mapstd(diab(:,1:8)')';
diabTrain = diab1(1:575,:)';
diabTest = diab1(576:768,:)';
diabTrainL = diab(1:575,9);
diabTestL = diab(576:768,9);
failscript
loadMNISTImages('train-images-idx3-ubyte');
trainImage = ans;
trainLabel = loadMNISTLabels('train-labels-idx1-ubyte');
testImage = loadMNISTImages('t10k-images-idx3-ubyte');
testLabel = loadMNISTLabels('t10k-labels-idx1-ubyte');
failscript
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier(trainImage, trainLabel, 1, 60000, 1);
failscript
waveform = waveform(randperm(size(waveform,1)),:);
waveformTrain = waveform(1:3750, 1:40)';
waveformTrainL = waveform(1:3750, 41);
waveformTest = waveform(3751:5000, 1:40)';
waveformTestL = waveform(3751:5000, 41);
failscript
%-- 14.04.2014 12:44 --%
load('local_uniform_2d.mat')
labels = zero(1,240000);
labels = zeros(1,240000);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier( Data, labels, 600, 400, 0.5);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier( Data, labels, 600, 400, 2);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier( Data, labels, 600, 400, 1);
[ nodes, edges, node_classes ] = construct_fastEGNG_classifier( Data, labels, 600, 400, 0.8);
%-- 14.04.2014 17:41 --%
ion = ion(randperm(size(ion,1)),:);
ionTrain = ion(1:260,1:34)';
ionTrainL = ion(1:260,35);
ionTestL = ion(261:351,35);
ionTest = ion(261:351,1:34)';
failscript
%-- 14.04.2014 23:20 --%
ion = ion(randperm(size(ion,1)),:);
ionTrain = ion(1:260,1:34)';
ionTrainL = ion(1:260,35);
ionTestL = ion(261:351,35);
ionTest = ion(261:351,1:34)';
failscript
cancer = cancer(randperm(size(cancer,1)),:);
cancerTrain = cancer(1:500,1:9)';
cancerTrainL = cancer(1:500,10);
cancerTest = cancer(501:683,1:9)';
cancerTestL = cancer(501:683,10);
failscript
newthyroid = newthyroid(randperm(size(newthyroid,1)),:);
newthyroid1 = mapstd(newthyroid(:,2:6)')';
thyrTrain = newthyroid1(1:160,:)';
thyrTest = newthyroid1(161:215,:)';
thyrTrainL = newthyroid(1:160,1);
thyrTestL = newthyroid(161:215,1);
failscript
