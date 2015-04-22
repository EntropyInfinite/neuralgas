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

Dims = size(Data, 2);

close all;
colordef white
pause on

% Step.0 Start with the root
agnh.NumOfNodes = 2;
agnh.ActualNodes = 2;
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
 deletionList = [];
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
        agnh.ActualNodes = agnh.ActualNodes+1;
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
     if agnh.nodes(k).deleted
         continue
     end
     %if there are no buffers, shrink the lambda or fix the node if
     %shrinking not necessary
     if isempty(agnh.nodes(k).buffer)
         if ~isempty(agnh.nodes(k).next_layer)
             agnh.nodes(k).fixed = true;
             continue
         end
         if ~agnh.nodes(k).fixed && isempty(agnh.nodes(k).maxdist)
            deletionList = [deletionList k];
            agnh.nodes(k).deleted = true;
            continue 
         end
         if (agnh.nodes(k).lambda>agnh.nodes(k).maxdist)
            agnh.nodes(k).lambda = agnh.nodes(k).maxdist;
         else
            agnh.nodes(k).fixed = true;
         end
     %else iterate over buffers and make new nodes out of them    
     else
        for l=1:length(agnh.nodes(k).buffer)
            agnh.TotalNodes = agnh.TotalNodes+1;
            agnh.ActualNodes = agnh.ActualNodes+1;
            agnh.nodes = [agnh.nodes AGNH_NewNode(agnh.nodes(k).buffer(l).mean, 2*max(agnh.nodes(k).buffer(l).std))];
            agnh.nodes(k).next_layer = [agnh.nodes(k).next_layer agnh.TotalNodes];
        end
        agnh.nodes(k).buffer = [];
     end
 end
 
 length(deletionList)
 
 for k = deletionList
     agnh.ActualNodes = agnh.ActualNodes-1;
     if any(agnh.next_layer==k)
         agnh.next_layer(agnh.next_layer==k) = [];
         continue
     end
     for i=1:length(agnh.nodes)
         if any(agnh.nodes(i).next_layer==k)
             agnh.nodes(i).next_layer(agnh.nodes(i).next_layer==k) = [];
             break
         end
     end
 end
 
 AGNH_GlobalPlot(agnh,Data)
 
 if all([agnh.nodes(~[agnh.nodes.deleted]).('fixed')])
     converged = true;
     iterCount
 end
 %TODO
 
 %shrink existing lambdas until possible
 %if not possible to shrink - fix it
 % if all fixed - end algorithm
end

return
end
