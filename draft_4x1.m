o = EFMouse; 
o.ID = '4x1';  
o.dir = '/Users/rubensanchez/desktop/EFMouse/new/4x1';
o.log  = true;  % Create a log file.
o.initialize(overwrite=true); 

%% add electrodes
o.addElectrode(tag = "Anterior",current = 0.05,type="surface",shape="circle",center= [-3.5,30,6],radius=0.6, thickness=1);
o.addElectrode(tag = "Posterior",current = 0.05,type="surface",shape="circle",center= [-3.5,25,6],radius=0.6, thickness=1);
o.addElectrode(tag = "Lateral",current = 0.05,type="surface",shape="circle",center= [-6,27.5,4],radius=0.6, thickness=1);
o.addElectrode(tag = "Medial",current = 0.05,type="surface",shape="circle",center= [-1,27.5,6],radius=0.6, thickness=1);
o.addElectrode(tag = "Center",current = -0.2,type="surface",shape="circle",center= [-3.5,27.5,6],radius=0.6, thickness=1);

o.run(targetStage=Stage.MESH,show=true);
ylim([-37 37])

%% export mesh
o.run(targetStage=Stage.EXPORT)

%% run FEM model

o.run(targetStage=Stage.GETDP,show=false)