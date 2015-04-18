function agnh = AGNH( Data, classLabels, NumOfSamples)
%AGNH structure creation: provides AGNH in a struct

% Check label data
if size(classLabels, 2) > 1
    if size(classLabels, 1) > 1
        error('Incorrect label data format!');
    else
        classLabels = classLabels';
    end
end

% Check data sample number
if NumOfSamples > size(Data, 1)
    error('Not enough data provided!');
end

close all;
colordef white
pause on

% Step.0 Start with the root
agnh.NumOfNodes = 2;
agnh.TotalNodes = 2;
agnh.nodes(1) = AGNH_NewNode(Data(1,:), sqrt(sum((Data(1,:)-Data(2,:)).^2)));
agnh.nodes(2) = AGNH_NewNode(Data(2,:), sqrt(sum((Data(1,:)-Data(2,:)).^2)));
% agnh.nodes(1).coord = Data(1,:);
% agnh.nodes(2).coord = Data(2,:);
% agnh.nodes(1).lambda = Inf;
% agnh.nodes(2).lambda = Inf;
agnh.next_layer = [1 2];
agnh.buffer = [];

n = 3;
converged = false;
iterCount = 0;

% epochs loop
while converged==false
 iterCount = iterCount+1;   
 % samples loop   
 for n=n:size(Data,1)
    curSample = Data(n,:);
    curNode = agnh;
    curIndex = -1;
    done = false;
    
    %find winner node
    while done == false
        if isempty(curNode.next_layer)
            nodes_containing = false;
        else
            sqdists = zeros(1,length(curNode.next_layer));
            i = 1;
            for k = curNode.next_layer
                sqdists(i) = sum((agnh.nodes(k).coord - curSample).^2);
                i = i+1;
            end
            nodes_containing = sqrt(sqdists)<=[agnh.nodes(curNode.next_layer).('lambda')];
        end
        
        if any(nodes_containing) == false
            %if no further winner node, algorithm proceeds to the next step
            done = true;
        elseif sum(nodes_containing)>1
            % if more than one winning node, choose closest one
            [~, I] = min(sqdists);
            curIndex = curNode.next_layer(I);
            curNode = agnh.nodes(curIndex);
        else
            % only one node wins, pick it
            curIndex = curNode.next_layer(nodes_containing);
            curNode = agnh.nodes(curIndex);
        end
    end
    
    if curIndex<0
        agnh.NumOfNodes = agnh.NumOfNodes+1;
        agnh.TotalNodes = agnh.TotalNodes+1;
        agnh.nodes = [agnh.nodes AGNH_NewNode(curSample, min(sqrt(sqdists)))];
        agnh.next_layer = [agnh.next_layer agnh.TotalNodes];
        continue;
    end
    
    curSampleDist = sqrt(sum((curSample-curNode.coord).^2));
    %if size(Data,2) == 2 || size(Data,2) == 3
    %    AGNH_Plot(agnh, Data, curIndex, curSample);
    %end
    
    % if there are no buffers
    if isempty(curNode.buffer)
        % if sample is farther than half lambda
        % assign it to the new buffer
        if curSampleDist>curNode.lambda/2
            tmp.points = curSample;
            tmp.mean = curSample;
            tmp.std = 0;
            curNode.buffer = [curNode.buffer tmp];
        % else update the actual used area if necessary
        else
            if isempty(curNode.maxdist) || (curNode.maxdist<curSampleDist)
                curNode.maxdist = curSampleDist;
            end
        end
    else
    % check the distance to all the buffers (including the node itself)
        buffDists = zeros(1,length(curNode.buffer));
        for k = 1:length(curNode.buffer)
            buffDists(k) = sqrt(sum((curSample-curNode.buffer(k).mean).^2));
        end
        if curSampleDist>curNode.lambda/2 && all(buffDists>curNode.lambda/2)
            tmp.points = curSample;
            tmp.mean = curSample;
            tmp.std = 0;
            curNode.buffer = [curNode.buffer tmp];
            buffIndex = -1;
        elseif any(buffDists<curSampleDist)
            [~, buffIndex] = min(buffDists);
        else
            buffIndex = -1;
        end
        if buffIndex>0
            curNode.buffer(buffIndex).points = [curNode.buffer(buffIndex).points; curSample];
            curNode.buffer(buffIndex).mean = mean(curNode.buffer(buffIndex).points);
            curNode.buffer(buffIndex).std = std(curNode.buffer(buffIndex).points);
        end
    end
    agnh.nodes(curIndex) = curNode;
 end
 % reset sample counter
 n = 1;
 %do stuff at the end of epoch
 
 %iterate through all the nodes
 for k=1:length(agnh.nodes)
     %if there are no buffers, shrink the lambda or fix the node if
     %shrinking not necessary
     if isempty(agnh.nodes(k).buffer)
         if ~isempty(agnh.nodes(k).next_layer)
             agnh.nodes(k).fixed = true;
             continue
         end
         if ~isempty(agnh.nodes(k).maxdist) && (agnh.nodes(k).lambda>agnh.nodes(k).maxdist)
            agnh.nodes(k).lambda = agnh.nodes(k).maxdist;
         else
            agnh.nodes(k).fixed = true;
         end
     %else iterate over buffers and make new nodes out of them    
     else
        for l=1:length(agnh.nodes(k).buffer)
            agnh.TotalNodes = agnh.TotalNodes+1;
            agnh.nodes = [agnh.nodes AGNH_NewNode(agnh.nodes(k).buffer(l).mean, 2*max(agnh.nodes(k).buffer(l).std))];
            agnh.nodes(k).next_layer = [agnh.nodes(k).next_layer agnh.TotalNodes];
        end
        agnh.nodes(k).buffer = [];
     end
 end
 
 fixedCount = sum([agnh.nodes.('fixed')])
 agnh.TotalNodes
 currentRatio = sum([agnh.nodes.('fixed')])/agnh.TotalNodes
 if all([agnh.nodes.('fixed')])
     converged = true;
     iterCount
 end
 %TODO
 
 %shrink existing lambdas until possible
 %if not possible to shrink - fix it
 % if all fixed - end algorithm
end

return
%ClassIdents = unique(classLabels);node
%NumOfClasses = numel(ClassIdents);
% SecondNodeIndex = find(classLabels~=classLabels(1,1),1);
% nodes = [Data(:,1)  Data(:,SecondNodeIndex)];
fgng.nodes = [Data(1,:)  Data(2,:)];
%classfreqs = [zeros(NumOfClasses,1) zeros(NumOfClasses,1)];
lambda = 2;
%lambda = norm(fgng.nodes(1)-fgng.nodes(2))/2;
fgng.node_lambdas = [lambda lambda];
fgng.point_coverages = [1 1];
fgng.node_classes = [-1 -1];
fgng.fixed_nodes = uint8([0 0]);
failed = 0;

% Initial connections (edges) matrix.
fgng.edges = [0  1;
                1  0;];
     
% Initial ages matrix.
ages = [ NaN  0;
               0  NaN;];

% Initial err Vector.
% err = [0 0];

Cur_NumOfNodes = 0;
Epoch = 0;

% scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/2 scrsz(4)/3-50 scrsz(3)/2 2*scrsz(4)/3])
ecount = 0;
plotArray = [2; 0]; %prev point for node count graph
tic
n = 3;
while ~all(fgng.fixed_nodes)
    ecount = ecount+1;
%     if toc>90
%         %dbstop in construct_fastEGNG_classifier2 at 93;
%         tic;
%         failed = 1;
%         return
%     end
for kk=1:NumOfEpochs
    fgng.point_coverages = ones(1,NumOfNodes);
    
    % Choose the next Input Training Vectors.
    nextblock = (kk-1)*NumOfSamples+1:1:kk*NumOfSamples;
    In = Data(:,nextblock);
    
    %set initial lambda

    for n=n:NumOfSamples

    % Step.1 Generate an input signal � according to P(�).
    params(9) = n;
    Input = In(:,n);

    %% Step 2. Find the two nearest units s1 and s2 to the new data sample.
    [s1, s2] = findTwoNearest(Input,fgng.nodes);
    params(10) = s1;
    params(11) = s2;

    %% Steps 3-6. Increment the age of all edges emanating from s1 .
    
        % Step 3. Increment the age of all edges emanating from s1. 
    s1_Neighbors = find(fgng.edges(:,s1)==1);
    SizeOfNeighborhood = length(s1_Neighbors);

    ages(s1_Neighbors,s1) = ages(s1_Neighbors,s1) + age_inc;
    ages(s1,s1_Neighbors) = ages(s1_Neighbors,s1);

    % Step 4. Add the squared distance to a local error counter variable:
    %error(s1) = error(s1) + distances(s1)^2;

    % Step 5. Move s1 and its topological neighbors towards �.
    ndist = norm(fgng.nodes(:,s1)-Input);
    if fgng.fixed_nodes(1,s1) == 1
        if ndist > fgng.node_lambdas(1, s1)
                 % Add the new node at target input point: 
               fgng.nodes = [fgng.nodes Input];

               NumOfNodes = NumOfNodes+1;
               r = NumOfNodes;

               % Insert edges connecting the new unit r with units q anf f. 
               fgng.edges = [fgng.edges  zeros(NumOfNodes-1,1)];
               fgng.edges = [fgng.edges; zeros(1,NumOfNodes)];
               fgng.edges(s1,r) = 1;
               fgng.edges(r,s1) = 1;
               ages = [ages  NaN*ones(NumOfNodes-1,1)];
               ages = [ages; NaN*ones(1,NumOfNodes)];
               ages(s1,r) = 0;
               ages(r,s1) = 0;
               fgng.fixed_nodes = [fgng.fixed_nodes 0];
               fgng.node_classes = [fgng.node_classes -1];
               fgng.node_lambdas = [fgng.node_lambdas ndist/2];
               fgng.point_coverages = [fgng.point_coverages 1];
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
                fgng.edges(s1,s2) = 1;
                fgng.edges(s2,s1) = 1;
                ages(s1,s2) = 0;
                ages(s2,s1) = 0;
                fgng.point_coverages(1,s1) = fgng.point_coverages(1,s1)+1;
            %end
        end
    else
        fgng.nodes(:,s1) = fgng.nodes(:,s1) + eb*(Input-fgng.nodes(:,s1));
        fgng.nodes(:,s1_Neighbors) = fgng.nodes(:,s1_Neighbors) + en*(repmat(Input,[1 SizeOfNeighborhood])-fgng.nodes(:,s1_Neighbors));

        fgng.edges(s1,s2) = 1;
        fgng.edges(s2,s1) = 1;
        ages(s1,s2) = 0;
        ages(s2,s1) = 0;

        if norm(fgng.nodes(:,s1)-Input) <= fgng.node_lambdas(1, s1)
            fgng.fixed_nodes(1,s1) = 1;
            fgng.point_coverages(1,s1) = fgng.point_coverages(1,s1)+1;
            fgng.node_classes(1,s1) = classLabels(1,n);
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
    [DelRow, DelCol] = find(ages>max_age);
    SizeDeletion = length(DelRow);
    for i=1:SizeDeletion
        fgng.edges(DelRow(i),DelCol(i)) = 0;
          ages(DelRow(i),DelCol(i)) = NaN;
    end
    %[nodes, edges, ages, fixed_nodes, node_classes] = edgeManagementF(In(:,n),nodes,edges,ages,fixed_nodes,node_classes,distances, params);

    %% Step 7. Dead Node Removal Procedure. 
    ii = 1;
    while NumOfNodes >= ii
        if any(fgng.edges(ii,:)) == 0

            fgng.edges = [fgng.edges(1:ii-1,:); fgng.edges(ii+1:NumOfNodes,:);];
            fgng.edges = [fgng.edges(:,1:ii-1)  fgng.edges(:,ii+1:NumOfNodes);];

            ages = [ages(1:ii-1,:); ages(ii+1:NumOfNodes,:);];
            ages = [ages(:,1:ii-1)  ages(:,ii+1:NumOfNodes);];

            fgng.nodes = [fgng.nodes(:,1:ii-1) fgng.nodes(:,ii+1:NumOfNodes);];
            fgng.fixed_nodes = [fgng.fixed_nodes(1,1:ii-1) fgng.fixed_nodes(1,ii+1:NumOfNodes);];
            fgng.node_classes = [fgng.node_classes(1,1:ii-1) fgng.node_classes(1,ii+1:NumOfNodes);];
            fgng.node_lambdas = [fgng.node_lambdas(1,1:ii-1) fgng.node_lambdas(1,ii+1:NumOfNodes);];
            fgng.point_coverages = [fgng.point_coverages(1,1:ii-1) fgng.point_coverages(1,ii+1:NumOfNodes);];

            NumOfNodes = NumOfNodes - 1;

            ii = ii -1; 
        end
        ii = ii+1;
    end
    %[nodes, edges, ages, fixed_nodes, node_classes] = removeUnconnectedF(nodes, edges,ages,fixed_nodes,node_classes);
    %%
    end
    n = 1;
    %hold on;
    plotArray = [plotArray [NumOfNodes; sum(fgng.fixed_nodes>0)]];
    %plot(0:ecount, plotArray(1,:), 0:ecount, plotArray(2,:));
    %plot(0:ecount, 1-plotArray(2,:)./plotArray(1,:));
    drawnow;
    

     %% Refresh drawing buffers
%     NumOfNodes = size(nodes,2);
    Cur_NumOfNodes = [Cur_NumOfNodes NumOfNodes];
    if length(Cur_NumOfNodes)>100
        Cur_NumOfNodes = Cur_NumOfNodes(end-100:end);
    end
% 
%     Cur_RMSE(kk) = norm(err)/sqrt(NumOfNodes);
%     RMSE = [RMSE Cur_RMSE(kk)];
%     if length(RMSE)>100
%         RMSE = RMSE(end-100:end);
%     end
% 
    Epoch = [Epoch ecount];
    if length(Epoch)>100
        Epoch = Epoch(end-100:end);
    end

      %% Plot everything
     subplot(1,2,1);
     hold on
     scatter(Data(1,:),Data(2,:),8,'y');
     %scatter3(Data(1,:),Data(2,:),Data(3,:),8,'y');
     plotgng2015(fgng.nodes,fgng.edges,'n',fgng.fixed_nodes);
     %viscircles(fgng.nodes',fgng.node_lambdas,'EdgeColor','r','LineWidth',1);
     hold off
     % xlim([-1/2 2.5]);
     % ylim([-1 8]);
     % zlim([-1/2 1.5]);
     % xlim([-1 6]);
     % ylim([-1 6]);
     % zlim([-7 7]);
     drawnow;
%  
      subplot(2,2,2);
      plot(1:ecount, 1-plotArray(2,2:end)./plotArray(1,2:end));
%      plot(Epoch,0,'r.');
      title('Free nodes proportion');
%      if kk>100
%           xlim([Epoch(1) Epoch(end)]);
%      end
      xlabel('Training Epoch Number');
      grid on;
%  
     subplot(2,2,4);
     plot(Epoch,Cur_NumOfNodes,'g.');
     title('Number of Neural Units in the Growing Neural Gas');
     if kk>100
       xlim([Epoch(1) Epoch(end)]);
     end
     xlabel('Training Epoch Number');
     grid on;
     M(kk)=getframe(gcf);
end

    if size(fgng.nodes,2)>= NumOfEpochs*NumOfSamples
        disp(lambda);
        error('Neural gas did not converge!');
    end
end
%NumOfNodes = size(nodes,2);


%node_classes = classLabels(:,nearest_indices);


end
