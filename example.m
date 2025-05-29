%% Compact example. 
o = EFMouse; % Create a default empty object of the class EFMouse
o.ID = 'retest';  % A name/tag for this simulation.
o.dir = '/Users/rubensanchez/desktop/EFMouse/1x1Montage_retest';
o.log  = true;  % Create a log file.
o.initialize(overwrite=true); 
% Define electrodes and craniotomy on the posterior left hemishphere
% Use surface electrodes 
o.addElectrode(tag = "Anterior",current = 0.2,type="surface",shape="circle",center= [-3.5,30,6],radius=0.6, thickness=1);
o.addElectrode(tag = "Posterior",current = -0.2,type="surface",shape="circle",center= [-3.5,24,6],radius=0.6, thickness=1);
o.addCraniotomy(tag = "left",center=[-3.5,27,6],radius=1.5,material=["csf","csf"]);
% Run meshing and simulation 


%o.run(targetStage=Stage.GETDP,show=false);

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
