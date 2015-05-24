%% AGNH plot
function AGNH_GlobalPlot(agnh, Data, classLabels)
    hold on
    scatter(Data(:,1),Data(:,2),5,'g');
    node_coords = reshape([agnh.nodes(~[agnh.nodes.deleted]).coord],size(Data,2),[])';
    scatter(node_coords(:,1),node_coords(:,2),6,'b');
    node_lambdas = reshape([agnh.nodes(~[agnh.nodes.deleted]).lambda],[],1);
    %if size(node_coords,2)>2
    %    size(node_coords)
    %end
    viscircles(node_coords(:,1:2),node_lambdas,'LineWidth',1);
    fixedCount = sum([agnh.nodes(~[agnh.nodes.deleted]).('fixed')]);
    currentRatio = sum([agnh.nodes(~[agnh.nodes.deleted]).('fixed')])/agnh.ActualNodes;
    titlestr = ['Fixed nodes ' num2str(fixedCount) '/' num2str(agnh.ActualNodes) ' = ' num2str(currentRatio)];
    title(titlestr)
    %pause
    hold off
end