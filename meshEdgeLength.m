% can adapt this for any mesh
mouse = EFMouse(dir='/Users/rubensanchez/desktop/EFMouse/4x1Montage_rubber',ID='4x1_rubber');
node = mouse.mesh.node';
elem = mouse.mesh.elem';

% each row in elem is 1 tetrahedron defined by 4 nodes
% a tetrahedron has 6 edges
% edge_all is a 2-D matrix encoding edges as pair of nodes (node_1,node_2) 
% for all tetrahedra in elem. 
edges_all = [
    elem(:,[1 2]);
    elem(:,[1 3]);
    elem(:,[1 4]);
    elem(:,[2 3]);
    elem(:,[2 4]);
    elem(:,[3 4])];
% make sure that (i,j) and (j,i) both are expressed as (i,j)
edges_sorted = sort(edges_all, 2);
% and then to not repeat edges
edges = unique(edges_sorted, 'rows');

% Get node spatial coordinates defined in node
% each row in node are x,y,z coordinates of a node
p1 = node(edges(:,1), :);  % Start point of edge
p2 = node(edges(:,2), :);  % End point of edge

%Euclidean distance
edge_lengths = sqrt(sum((p2 - p1).^2, 2));

% to check the shape of the distribution
histogram(edge_lengths);
min_edge_length = min(edge_lengths);
max_edge_length = max(edge_lengths);
% use median length for the rubber electrodes.
median_edge_length = median(edge_lengths);
disp(median_edge_length)
