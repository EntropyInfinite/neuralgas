function predicted = AGNHc_classify(agnh, newSet)
    
    predicted = zeros(size(newSet,1),length(agnh.Classes));
    cc = 0;
    
    for m=1:size(newSet,1)
        
        probcounter = zeros(1,length(agnh.Classes));        
        curSample = newSet(m,:);
        curNode = agnh;
        curIndex = -1;
        
        done = false;
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
                if curIndex<0
                    [~, I] = min(sqdists);
                    curIndex = curNode.next_layer(I);
                    curNode = agnh.nodes(curIndex);
                    probcounter = probcounter + curNode.classcount;
                    cc = cc+1;
                end
                done = true;
            elseif sum(nodes_containing)>1
                % if more than one winning node, choose closest one
                [~, I] = min(sqdists);
                curIndex = curNode.next_layer(I);
                curNode = agnh.nodes(curIndex);
                probcounter = probcounter + curNode.classcount;
            else
                % only one node wins, pick it
                curIndex = curNode.next_layer(nodes_containing);
                curNode = agnh.nodes(curIndex);
                probcounter = probcounter + curNode.classcount;
            end
        end
        
        probcounter = probcounter/sum(probcounter);
        predicted(m,:) = probcounter;
        
    end
    cc
end

