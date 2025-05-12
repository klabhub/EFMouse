% right now no gel layers. Need to be added

% run this after running test_for_rubber_electrode.m
o = EFMouse(dir='/Users/rubensanchez/desktop/EFMouse/4x1Montage_border',ID='4x1_border');

mouse_node = o.mesh.node';
mouse_elem = o.mesh.elem';
mouse_label = o.mesh.label';

% make a findElement electrode and then copy it and move it up

i = 5; % the other four are stimulation electrodes already defined
o.eCenter(:,i) = [-1.1819,-9.655,10.129]; %this is just aux, not the final center coordinates of the electrode
o.eRadius(:,i) = 12.25/2; 
aux_electrode = findElements(o.model.Mesh,'radius',...
                        o.eCenter(:,i),...
                        o.eRadius(i));

tissue_touched = unique(o.mesh.label(aux_electrode));
elem_tiss_touched  = zeros(1,numel(tissue_touched));

fprintf(' aux electrode: %s: touching tissue:\n',o.eTag(i));
for t = 1:numel(tissue_touched)
    tiss = tissue_touched(t);
    elem_tiss_touched(1,t) = sum(o.mesh.label(aux_electrode) == tiss);
    fprintf('   %s: num elements = %d\n', o.labelToTissue(tiss),elem_tiss_touched(1,t));
end

% this is to define the electrode only for the skin elements
[~,idx] = max(elem_tiss_touched);
tiss_elec = tissue_touched(idx);
fprintf('   Creating %s aux electrode %s in max touched tissue: %s.\n',o.eShape(i),o.eTag(i),o.labelToTissue(tiss_elec))
aux_electrode = aux_electrode(o.mesh.label(aux_electrode)==tiss_elec);
aux_elem = mouse_elem(aux_electrode,:);


% create the "floating" electrode
% This was the first try. It did not work, because we only want a "thin"
% rubber electrode.

% find all the nodes ids
aux_node_idx = unique(aux_elem);
% find the coordinates for those nodes
aux_node = mouse_node(aux_node_idx,:);

% translate the nodes in the z direction (+z) ~pz mm
pz = 6;
elec_node = aux_node;
elec_node(:,3) = aux_node(:,3) + pz;

% get number of nodes for re-indexing
num_mouse_node = size(mouse_node, 1);
% append new nodes to the mouse nodes 
comb_node = [mouse_node; elec_node];

% create node IDs starting from the last node ID of the original mouse mesh
% eg. for N nodes, N+1, N+2, etc.
elec_node_idx = [num_mouse_node + (1:length(aux_node_idx))]';
% create a dictionary to map from old node id to new node id (the ones of
% the final electrode)
node_mapping = dictionary(aux_node_idx, elec_node_idx);

% update node ids using the node_mapping dictionary
elec_elem = node_mapping(aux_elem); 

% append new elements
comb_elem = [mouse_elem; elec_elem];

% add elem labels
elec_label = ones(size(elec_elem,1),1)*12;

comb_label = [mouse_label;elec_label];

% redefine the mesh using the added electrode
o.mesh.node = comb_node';
o.mesh.elem = comb_elem';
o.mesh.label = comb_label';


% plot and save as necessary
%pdemesh(comb_node',comb_elem',FaceColor='red');
%save(file(o,"OBJECT"),"o");