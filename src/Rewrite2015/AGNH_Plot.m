%% AGNH plot
function AGNH_Plot(agnh, Data, curIndex, curSample)
    scatter(Data(:,1),Data(:,2),2,'g');
    hold on
    node_coords = reshape([agnh.nodes.coord],2,[])';
    scatter(node_coords(:,1),node_coords(:,2),6,'b');
    if curIndex>0
        scatter(node_coords(curIndex,1),node_coords(curIndex,2),12,'r');
    end
    scatter(curSample(:,1),curSample(:,2),12,'k');
    node_lambdas = reshape([agnh.nodes.lambda],[],1);
    if size(node_coords,2)>2
        size(node_coords)
    end
    viscircles(node_coords,node_lambdas,'LineWidth',1);
    pause(0.2)
    hold off
end