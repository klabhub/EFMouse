% Tetrahedral elements (4 nodes each)
faces = [elements(:,[1 2 3]);
         elements(:,[1 2 4]);
         elements(:,[1 3 4]);
         elements(:,[2 3 4])];
 
% Sort faces to ignore orientation
sortedFaces = sort(faces, 2);
 
% Count face appearances
[uniqueFaces, ~, ic] = unique(sortedFaces, 'rows');
counts = accumarray(ic, 1);
 
% Boundary faces appear only once
boundaryFaces = uniqueFaces(counts == 1, :);
 
% Boundary nodes
boundaryNodes = unique(boundaryFaces(:));