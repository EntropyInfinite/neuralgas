function newNode = AGNH_NewNode( nodeCoords, lambda )
    newNode.coord = nodeCoords;
    newNode.lambda = lambda;
    newNode.next_layer = [];
    newNode.buffer = [];
    newNode.maxdist = [];
    newNode.fixed = false;
    newNode.deleted = false;
    
end

