function  [nodes,edges,ages,fixed_nodes,error] = addNewNeuronF(Input,s1,nodes,edges,ages,fixed_nodes,error)  %#codegen 
                                                                                                                                    % Checks whether this function is
                                                                                                                                    % suitable for automatic .mex code generation.
% Add the new node at target input point: 
   nodes = [nodes Input];
   
   NumOfNodes = size(nodes,2);
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
   error = [error -1];

% MEX-code Generation:  
   
% mexcfg = coder.config('mex');
% mexcfg.DynamicMemoryAllocation = 'AllVariableSizeArrays'; 
% codegen -config mexcfg addNewNeuron.m -args {coder.typeof(nodes,[Inf Inf]), coder.typeof(edges,[Inf Inf]), coder.typeof(ages,[Inf Inf]), coder.typeof(error,[1 Inf]), double(0)}
