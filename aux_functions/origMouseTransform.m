% function to download, reorient and relabel the Alekseichuk et al. mouse 
% reference data website: https://zenodo.org/record/3857041#.ZBUKEXaZOUn
% and paper: https://www.sciencedirect.com/science/article/pii/S1053811919302320
% which itself was modified from the original Digimouse
% Digimouse: https://neuroimage.usc.edu/neuro/Digimouse

% the *.mat file was build by us from the original Alekseichuk download *.msh file    

load('aux_files/Mouse_Digimouse.mat')
node = msh.POS'; % nodes
elem = msh.TETS'; % tetrahedra
face = msh.TRIANGLES'; % faces

%% Relabel tissue.This is to get tissue labels in order 1,2,3,etc
% from the Alekseichuk et al., labels:
% Grey matter (2)-->(1)
% CSF (3)-->(2)
% Bone volume (4)-->(3)
% Soft tissues (5)-->(4)
% Eyeballs volume (8)-->(5)
old_elem_labels = [2,3,4,5,8];
new_elem_labels = [1,2,3,4,5];      
elem(5,:) = changem(elem(5,:),new_elem_labels,old_elem_labels);

%% Orientation correction. To match to stereotaxic orientation
% create a node_aux variable that will reflect any changes make to the
% original
% flip x node coordinates (left -->+ right: x axis)
node_aux(1,:) = -node(1,:);
% exchange y and z coordinates.
node_aux(2,:) = node(3,:); % (posterior -->+ anterior: y axis)
node_aux(3,:) = node(2,:); % (inferior -->+ superior: z axis)
% the modified node is the new node
node = node_aux;

%% Save the modified mesh
% This is the mesh used for our ef stimulation modeling
save('aux_files/EFMouse_mesh_clean.mat','node','elem','face');



