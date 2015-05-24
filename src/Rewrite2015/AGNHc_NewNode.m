function newNode = AGNHc_NewNode( nodeCoords, lambda, classcount )
    newNode.coord = nodeCoords;
    newNode.lambda = lambda;
    newNode.next_layer = [];
    newNode.buffer = [];
    newNode.maxdist = [];
    newNode.fixed = false;
    newNode.deleted = false;
    newNode.class = [];
    newNode.classcount = zeros(1,classcount);
end

