% test for merging electrodes
%% define stim electrodes and craneotomy
o = EFMouse; 
o.dir = '/Users/rubensanchez/desktop/EFMouse/4x1Montage_squareMerge'; 
o.ID = '4x1_squareMerge';  
o.log  = true; 
o.initialize(overwrite=true);

o.eTag = ["Anterior" "Posterior" "Lateral" "Medial" "Lumbar"];
o.eShape = ["circular" "circular" "circular" "circular" "circular"];
o.eCurrent = [0.05,0.05,0.05,0.05,-0.2];
o.eCenter = [-3.56,29.5,5.45;
             -3.43,24.8,6.02;
             -5.5,26.72,3.58; 
              -1.2,27.11,6.29]';
o.eRadius = [0.71,0.67,0.64,0.6]';
o.cCenter = [-3.4236,27.1067,5]';
o.cRadius = 1.5;
o.run(targetStage=1,show=true)

% just for easier reference
mouse_node = o.mesh.node';
mouse_elem = o.mesh.elem';
mouse_label = o.mesh.label';


%% create the square electrode

rubber_elec = createpde();
x_size = 24.5/2;
y_size = 24.5/2;
z_size = 1;
gm = multicuboid(x_size, y_size, z_size);
rubber_elec.Geometry = gm;
edge_length = 0.3;
generateMesh(rubber_elec, Hmax=edge_length, GeometricOrder="linear");
rubber_node = rubber_elec.Mesh.Nodes';        
rubber_elem = rubber_elec.Mesh.Elements'; 
% position rubber electrode in chosen coordinates
px = -1.64209;
py = -14.5699;
pz = 8; %inside
rubber_node = rubber_node + [px,py,pz];

%% merge
tol = 0.1;  % adjust based on mesh scale (what is this?)
% Use KDTree to match nodes
Mdl = KDTreeSearcher(mouse_node);
[idx, dist] = knnsearch(Mdl, rubber_node);

% Nodes in mesh2 that match mesh1 (shared)
sharedIdx2 = find(dist < tol);
sharedIdx1 = idx(sharedIdx2);

% Unique nodes in rubber_node (not shared)
nonSharedIdx2 = [setdiff(1:size(rubber_node,1), sharedIdx2)]';

% Keep original mesh1 nodes
mergedNodes = mouse_node;

% Map from old mesh2 node indices to new merged node indices
nodeIDMap = containers.Map('KeyType','int32','ValueType','int32');

% Add non-shared nodes from rubber_mesh
for i = 1:length(nonSharedIdx2)
    newNodeID = size(mergedNodes, 1) + 1;
    mergedNodes(newNodeID,:) = rubber_node(nonSharedIdx2(i),:);
    nodeIDMap(nonSharedIdx2(i)) = newNodeID;
end

% Map shared nodes to existing mesh1 indices
for i = 1:numel(sharedIdx2)
    nodeIDMap(sharedIdx2(i)) = sharedIdx1(i);
end

rubber_elem_remapped = zeros(size(rubber_elem));
for k = 1:numel(rubber_elem)
    oldNode = rubber_elem(k);
    rubber_elem_remapped(k) = nodeIDMap(oldNode);
end

mergedElems = [mouse_elem; rubber_elem_remapped];

rubber_label = ones(size(rubber_elem_remapped,1),1)*12;

mergedLabels = [mouse_label;rubber_label];

pdemesh(mergedNodes',mergedElems')


% redefine the mesh using the added rubber electrode
o.mesh.node = mergedNodes';
o.mesh.elem = mergedElems';
o.mesh.label = mergedLabels';

o.eTissue(5) = 4;
o.tissueLabel("Lumbar") = 12;

o.tissueMaterial("Lumbar") = "conductor";

