% get the elements that are in the boundary
% Tetrahedral elements (4 nodes each)
aux_faces = [aux_elem(:,[1 2 3]);
         aux_elem(:,[1 2 4]);
         aux_elem(:,[1 3 4]);
         aux_elem(:,[2 3 4])];
 
% Sort faces to ignore orientation
sortedFaces = sort(aux_faces, 2);
 
% Count face appearances
[uniqueFaces, ~, ic] = unique(sortedFaces, 'rows');
counts = accumarray(ic, 1);
 
% Boundary faces appear only once
boundaryFaces = uniqueFaces(counts == 1, :);
 
% Boundary nodes
boundaryNodeIDs = unique(boundaryFaces(:));
boundaryNodes = mouse_node(boundaryNodeIDs, :);


