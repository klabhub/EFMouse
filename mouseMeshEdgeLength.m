mouse = EFMouse(dir='/Users/rubensanchez/desktop/EFMouse/4x1Montage_rubber',ID='4x1_rubber');
node = mouse.mesh.node';
elem = mouse.mesh.elem';

edges_all = [
    elem(:,[1 2]);
    elem(:,[1 3]);
    elem(:,[1 4]);
    elem(:,[2 3]);
    elem(:,[2 4]);
    elem(:,[3 4])];
% use this to make sure if (i,j) and (j,i) both are expressed as (i,j)
edges_sorted = sort(edges_all, 2);
edges = unique(edges_sorted, 'rows');

% Get node coordinates
p1 = node(edges(:,1), :);  % Start point of edge
p2 = node(edges(:,2), :);  % End point of edge

%Euclidean distance
edge_lengths = sqrt(sum((p2 - p1).^2, 2));