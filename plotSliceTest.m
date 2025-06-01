% user defined
slice = 0;
ax = "x";
% TODO: some errors if the slice selected is outside the mesh bounds
map_slice_dir = dictionary(["x","y","z"],[1,2,3]);
row = map_slice_dir(ax);
slice_node_idx = find(abs(o.mesh.node(row,:) - slice) <= 0.5);
mask = all(ismember(o.mesh.elem,slice_node_idx),1);
slice_elem = o.mesh.elem(:,mask);
slice_label = o.mesh.label(mask);
figure;
pdeplot3D(o.mesh.node,slice_elem(:,slice_label==1),ColorMapData = data);

%xlim([-4,6])
%ylim([slice-1,slice+1])
%zlim([0,6])
if strcmp(ax,"y")
    view([0 0])
elseif strcmp(ax,"x")
    view([90 0])
else
    view([0 90])
end