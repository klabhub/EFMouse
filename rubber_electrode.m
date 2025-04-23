%% Create rubber electrode mesh
rubber_elec = createpde();
% Following Sanchez-Leon et al.2025 paper: 600 mm^2 ~ 24.5 mm x 24.5 mm
% https://elifesciences.org/reviewed-preprints/100941v2#s2
% 600 mm^2 seems too big, divide by two (for now)
% check in more papers

% ****User defined:
x_size = 24.5/2;
y_size = 24.5/2;
z_size = 1;

% create a Matlab geometry object, discrete geometry
% https://www.mathworks.com/help/pde/ug/multicuboid.html
% https://www.mathworks.com/help/pde/ug/pde.discretegeometry.html
gm = multicuboid(x_size, y_size, z_size);
rubber_elec.Geometry = gm;
% tetrahedral mesh for a 3-D geometry 
% saved in rubber_elec.Mesh
% see mouseMeshEdgeLength.m to get edge_length value
% use the resolution of the mouse mesh for consistent modelling
% (for previous circular electrodes and craniotomy, this was not used)
%edge_length = 0.2188;
% TODO: use only the info. from the tissue used: skin.
% Right now the function takes all the mouse mesh into account and
% the brain resolution is higher so it affects the median value.
% the following seems to approximate the skin mesh resolution.
edge_length = 0.3; % this seems to better approximate the skin
generateMesh(rubber_elec, Hmax=edge_length, GeometricOrder="linear"); %linear defines a 4 nodes tetrahedron. see generateMesh.m
% get nodes and elements
rubber_node = rubber_elec.Mesh.Nodes';        
rubber_elem = rubber_elec.Mesh.Elements';   

% assign label number for the rubber electrode.
% for now use the same label as circular electrodes.
% new info. needs to be added: rubber label and conductivity
% use label = 12 for lumbar return
% this is defined from number of tissue + number of electrodes
rubber_label = ones(size(rubber_elem,1),1)*12;

% plot for checking
figure;
pdemesh(rubber_elec)
% also can use pdeviz(rubber_electrode.Mesh)

%% Open mouse mesh and extract node and elements and labels
% open like this is just for the test
o = EFMouse(dir='/Users/rubensanchez/desktop/EFMouse/4x1Montage_rubber',ID='4x1_rubber');
mouse_node = o.mesh.node';
mouse_elem = o.mesh.elem';
mouse_label = o.mesh.label';
% plot 
figure;
pdemesh(mouse_node',mouse_elem',FaceColor='white');

% ****User defined:
% use the plot to get the xy,z center where the rubber electrode will be positioned
px = -1.64209;
py = -14.5699;
%pz = 8.92358;
%pz = 10; % some part touching
%pz = 11; % flying
pz = 8; %inside

%% move rubber electrode to the new position
% rubber electrode is created above in [0,0,0] center
rubber_node = rubber_node + [px,py,pz];
% plot for checking positioning
hold on;
pdemesh(rubber_node',rubber_elem')

%% combine the two meshes and reindex nodes
% concatenate the node coordinates info.
comb_node = [mouse_node; rubber_node];
% get number of nodes for reindexing
num_mouse_node = size(mouse_node, 1);
% change element node indices (each tetrahedron (row of elem) is defined 
% by 4 nodes), then concatenate elements.
comb_elem = [mouse_elem; (rubber_elem + num_mouse_node)];
% concatenate labels (they define the tissue and conductivity);
comb_label = [mouse_label;rubber_label];

% plot
figure;
pdemesh(comb_node',comb_elem',FaceColor='red');


% plot only label
figure;
pdemesh(comb_node',comb_elem(comb_label == 12,:)',FaceColor='green')


% redefine the mesh using the added rubber electrode
o.mesh.node = comb_node';
o.mesh.elem = comb_elem';
o.mesh.label = comb_label';

% add tissue and tissue label
o.eTissue(5) = 4;
o.tissueLabel("Lumbar") = 12;
% need to temporary change EFMouse to make this propery public, so it can
% be changed by hand. Check if error shows.
o.tissueMaterial("Lumbar") = "conductor";
