%% Create rubber electrode mesh
rubber_electrode = createpde();
% Followinf Sanchez-Leon et al.2025 paper: 600 mm^2 ~ 24.5 mm x 24.5 mm
% https://elifesciences.org/reviewed-preprints/100941v2#s2
x_size = 24.5;
y_size = 24.5;
z_size = 1;
% creates a Matlab geometry object, discrete geometry
% https://www.mathworks.com/help/pde/ug/multicuboid.html
% https://www.mathworks.com/help/pde/ug/pde.discretegeometry.html
gm = multicuboid(x_size, y_size, z_size);
rubber_electrode.Geometry = gm;
% linear tetrahedral mesh for a 3-D geometry (4 nodes x tetrahedron)
% saved in rubber_electrode.Mesh
generateMesh(rubber_electrode,Hmax=1, GeometricOrder="linear");
% can plot for checking
pdemesh(rubber_electrode)

%%