function [ nodes, edges, node_classes, node_lambdas, point_coverages, failed ] = adaptiveFEGNG_construct( Data, classLabels, NumOfEpochs, NumOfSamples)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 % Unsupervised Self Organizing Map. Growing Neural Gas (GNG) Algorithm.

% Main paper used for development of this neural network was:
% Fritzke B. "A Growing Neural Gas Network Learns Topologies", in 
%                         Advances in Neural Information Processing Systems, MIT Press, Cambridge MA, 1995.

% Check label data
if size(classLabels, 1) > 1
    if size(classLabels, 2) > 1
        error('Incorrect label data format!');
    else
        classLabels = classLabels';
    end
end

% Check data sample number
if NumOfEpochs*NumOfSamples > size(Data, 2)
    error('Not enough data provided!');
end

%clc;
%clear;
close all;
colordef black
%NumOfEpochs   = 600;
%NumOfSamples = 400;
age_inc               = 1;
max_age             = 100;
max_nodes         = 300;
eb                         = .05;
en                         = .005;
lambda                   = -0.5; %temporary value
alpha                    = .5;     % q and f units err reduction constant.
d                           = .99;   % err reduction factor.
RMSE                  = zeros(1,NumOfEpochs);

% load local_circular_2d1.mat
% load local_circular_2d2.mat;
% load local_sphere_shell.mat;
% load local_torus.mat;
% load local_uniform_2d.mat
% load local_quartersphere.mat;
% load multi_manifold_3d.mat;
 
% Define the params vector where the GNG algorithm parameters are stored:
params = [ age_inc;
                    max_age;
                    max_nodes;
                    eb;
                    en;
                    alpha;
                    lambda;
                    d;
                    0;   % Here, we insert the sample counter.
                    1;   % This is reserved for s1 (bmu);
                    2;]; % This is reserved for s2 (secbmu);

Cur_RMSE = zeros(1,NumOfEpochs);
RMSE = [];
Epoch = [];
Cur_NumOfNodes = [];

% Step.0 Start with two neural units (nodes) selected from input data:
NumOfNodes = 2;
ClassIdents = unique(classLabels);
NumOfClasses = numel(ClassIdents);
% SecondNodeIndex = find(classLabels~=classLabels(1,1),1);
% nodes = [Data(:,1)  Data(:,SecondNodeIndex)];
nodes = [Data(:,1)  Data(:,2)];
%classfreqs = [zeros(NumOfClasses,1) zeros(NumOfClasses,1)];
lambda = norm(nodes(1)-nodes(2))/2;
node_lambdas = [lambda lambda];
point_coverages = [1 1];
node_classes = [-1 -1];
fixed_nodes = uint8([0 0]);
failed = 0;

% Initial connections (edges) matrix.
edges = [0  1;
                1  0;];
     
% Initial ages matrix.
ages = [ NaN  0;
               0  NaN;];

% Initial err Vector.
% err = [0 0];

% scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/2 scrsz(4)/3-50 scrsz(3)/2 2*scrsz(4)/3])
tic
while ~all(fixed_nodes)
    if toc>60
        %dbstop in construct_fastEGNG_classifier2 at 93;
        tic;
        failed = 1;
        return
    end
for kk=1:NumOfEpochs
    point_coverages = ones(1,NumOfNodes);
    
    % Choose the next Input Training Vectors.
    nextblock = (kk-1)*NumOfSamples+1:1:kk*NumOfSamples;
    In = Data(:,nextblock);
    
    %set initial lambda

    for n=3:NumOfSamples

    % Step.1 Generate an input signal î according to P(î).
    params(9) = n;
    Input = In(:,n);

    %% Step 2. Find the two nearest units s1 and s2 to the new data sample.
    [s1 s2 distances] = findTwoNearest(Input,nodes);
    params(10) = s1;
    params(11) = s2;

    %% Steps 3-6. Increment the age of all edges emanating from s1 .
    
        % Step 3. Increment the age of all edges emanating from s1. 
    s1_Neighbors = find(edges(:,s1)==1);
    SizeOfNeighborhood = length(s1_Neighbors);

    ages(s1_Neighbors,s1) = ages(s1_Neighbors,s1) + age_inc;
    ages(s1,s1_Neighbors) = ages(s1_Neighbors,s1);

    % Step 4. Add the squared distance to a local error counter variable:
    %error(s1) = error(s1) + distances(s1)^2;

    % Step 5. Move s1 and its topological neighbors towards î.
    ndist = norm(nodes(:,s1)-Input);
    if fixed_nodes(1,s1) == 1
        if ndist > node_lambdas(1, s1)
                 % Add the new node at target input point: 
               nodes = [nodes Input];

               NumOfNodes = NumOfNodes+1;
               r = NumOfNodes;

               % Insert edges connecting the new unit r with units q anf f. 
               edges = [edges  zeros(NumOfNodes-1,1)];
               edges = [edges; zeros(1,NumOfNodes)];
               edges(s1,r) = 1;
               edges(r,s1) = 1;
               ages = [ages  NaN*ones(NumOfNodes-1,1)];
               ages = [ages; NaN*ones(1,NumOfNodes)];
               ages(s1,r) = 0;
               ages(r,s1) = 0;
               fixed_nodes = [fixed_nodes 0];
               node_classes = [node_classes -1];
               node_lambdas = [node_lambdas ndist/2];
               point_coverages = [point_coverages 1];
        else
%             if classLabels(1,n)~=node_classes(1,s1)
%                 nodes = [nodes Input];
%                 NumOfNodes = NumOfNodes+1;
%                 r = NumOfNodes;
%                 edges = [edges  zeros(NumOfNodes-1,1)];
%                 edges = [edges; zeros(1,NumOfNodes)];
%                 edges(s1,r) = 1;
%                 edges(r,s1) = 1;
%                 ages = [ages  NaN*ones(NumOfNodes-1,1)];
%                 ages = [ages; NaN*ones(1,NumOfNodes)];
%                 ages(s1,r) = 0;
%                 ages(r,s1) = 0;
%                 fixed_nodes = [fixed_nodes 1];
%                 node_classes = [node_classes classLabels(1,n)];
%                 node_lambdas(1,s1) = ndist*point_coverages(1,s1)/(point_coverages(1,s1)+1);
%                 node_lambdas = [node_lambdas ndist-node_lambdas(1,s1)];
%                 point_coverages = [point_coverages 1];
%                 point_coverages(1,s1) = point_coverages(1,s1)-1;
            %else
                edges(s1,s2) = 1;
                edges(s2,s1) = 1;
                ages(s1,s2) = 0;
                ages(s2,s1) = 0;
                point_coverages(1,s1)=point_coverages(1,s1)+1;
            %end
        end
    else
        nodes(:,s1) = nodes(:,s1) + eb*(Input-nodes(:,s1));
        nodes(:,s1_Neighbors) = nodes(:,s1_Neighbors) + en*(repmat(Input,[1 SizeOfNeighborhood])-nodes(:,s1_Neighbors));

        edges(s1,s2) = 1;
        edges(s2,s1) = 1;
        ages(s1,s2) = 0;
        ages(s2,s1) = 0;

        if norm(nodes(:,s1)-Input) <= node_lambdas(1, s1)
            fixed_nodes(1,s1) = 1;
            point_coverages(1,s1)=point_coverages(1,s1)+1;
            node_classes(1,s1) = classLabels(1,n);
        end
    end

    % Step 6.
    % If s1 and s2 are connected by an edge, set the age of this edge to zero.
    % If such an edge does not exist, create it.
    %edges(s1,s2) = 1;
    %edges(s2,s1) = 1;
    %ages(s1,s2) = 0;
    %ages(s2,s1) = 0;

    % Step 7. Remove edges with an age>max_age.
    [DelRow DelCol] = find(ages>max_age);
    SizeDeletion = length(DelRow);
    for i=1:SizeDeletion
        edges(DelRow(i),DelCol(i)) = 0;
          ages(DelRow(i),DelCol(i)) = NaN;
    end
    %[nodes, edges, ages, fixed_nodes, node_classes] = edgeManagementF(In(:,n),nodes,edges,ages,fixed_nodes,node_classes,distances, params);

    %% Step 7. Dead Node Removal Procedure. 
    ii = 1;
    while NumOfNodes >= ii
        if any(edges(ii,:)) == 0

            edges = [edges(1:ii-1,:); edges(ii+1:NumOfNodes,:);];
            edges = [edges(:,1:ii-1)  edges(:,ii+1:NumOfNodes);];

            ages = [ages(1:ii-1,:); ages(ii+1:NumOfNodes,:);];
            ages = [ages(:,1:ii-1)  ages(:,ii+1:NumOfNodes);];

            nodes = [nodes(:,1:ii-1) nodes(:,ii+1:NumOfNodes);];
            fixed_nodes = [fixed_nodes(1,1:ii-1) fixed_nodes(1,ii+1:NumOfNodes);];
            node_classes = [node_classes(1,1:ii-1) node_classes(1,ii+1:NumOfNodes);];
            node_lambdas = [node_lambdas(1,1:ii-1) node_lambdas(1,ii+1:NumOfNodes);];
            point_coverages = [point_coverages(1,1:ii-1) point_coverages(1,ii+1:NumOfNodes);];

            NumOfNodes = NumOfNodes - 1;

            ii = ii -1; 
        end
        ii = ii+1;
    end
    %[nodes, edges, ages, fixed_nodes, node_classes] = removeUnconnectedF(nodes, edges,ages,fixed_nodes,node_classes);
    %%
    end

    

     %% Refresh drawing buffers
%     NumOfNodes = size(nodes,2);
%     Cur_NumOfNodes = [Cur_NumOfNodes NumOfNodes];
%     if length(Cur_NumOfNodes)>100
%         Cur_NumOfNodes = Cur_NumOfNodes(end-100:end);
%     end
% 
%     Cur_RMSE(kk) = norm(err)/sqrt(NumOfNodes);
%     RMSE = [RMSE Cur_RMSE(kk)];
%     if length(RMSE)>100
%         RMSE = RMSE(end-100:end);
%     end
% 
%     Epoch = [Epoch kk];
%     if length(Epoch)>100
%         Epoch = Epoch(end-100:end);
%     end

      %% Plot everything
%      subplot(1,2,1);
%      plotgng(nodes,edges,'n');
%      % xlim([-1/2 2.5]);
%      % ylim([-1 8]);
%      % zlim([-1/2 1.5]);
%      % xlim([-1 6]);
%      % ylim([-1 6]);
%      % zlim([-7 7]);
%      drawnow;
%  
%      subplot(2,2,2);
%      plot(Epoch,0,'r.');
%      title('RMS err');
%      if kk>100
%           xlim([Epoch(1) Epoch(end)]);
%      end
%      xlabel('Training Epoch Number');
%      grid on;
%  
%      subplot(2,2,4);
%      plot(Epoch,Cur_NumOfNodes,'g.');
%      title('Number of Neural Units in the Growing Neural Gas');
%      if kk>100
%        xlim([Epoch(1) Epoch(end)]);
%      end
%      xlabel('Training Epoch Number');
%      grid on;
%      M(kk)=getframe(gcf);
end

    if size(nodes,2)>= NumOfEpochs*NumOfSamples
        disp(lambda);
        error('Neural gas did not converge!');
    end
end
%NumOfNodes = size(nodes,2);


%node_classes = classLabels(:,nearest_indices);
end
