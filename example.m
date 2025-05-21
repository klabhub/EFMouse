%% Compact example. 
o = EFMouse; % Create a default empty object of the class EFMouse
o.ID = 'retest';  % A name/tag for this simulation.
o.log  = true;  % Create a log file.
o.dir = 'c:/temp/efmouse'; % Results and the object (retest.mat) will be saved here.
o.initialize(overwrite=true); 
% Define electrodes and craniotomy on the posterior left hemishphere
% Use surface electrodes 
o.addElectrode(tag = "Anterior",current = -1,center= [-3  29 5],radius=0.7, thickness=1,type="surface",shape="circle");
o.addElectrode(tag = "Posterior",current = 1,center= [-3  25 5],radius=0.7, thickness=1, type="surface",shape="circle");
o.addCraniotomy(tag = "left",center=[-3 27 5],radius=1);
% Run meshing and simulation 
o.run(targetStage=Stage.GETDP,show=false);

%% Plot results

plotMesh(o) % The mesh with electrodes
plotEf(o,type='eMag',percentile=98,tissue='gray'); % Efield magnitude in gray matter
% Show field estimates for two tissue type that are part of the craniotomy
analyzeTissue(o,["leftskin" "leftbone"]);
% SHow field esimate for circular ROI
roi.shape = 'radius';
roi.center = [3 27 5];
roi.radius = 1;
analyzeRoi(o,roi);
