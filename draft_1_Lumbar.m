% notebook for the 1+lumbar montage.
% The stimulation electrode is positioned in the skin above visual cortex.
% The return electrode is a square positioned in the skin in the lumbar.

o = EFMouse; 
o.ID = '1xlumbar';  
o.dir = '/Users/rubensanchez/desktop/EFMouse/1xlumbar';
o.log  = true; 
o.initialize(overwrite=true); 

%% add electrodes
% for reference the center is located where the craniotomy eventually will be perfomed
o.addElectrode(tag = "Anterior",current = 0.2,type="surface",shape="circle", ...
    center= [-3.5,30,6],radius=0.6, thickness=1);
% the return electorde goes in the lumbar, test different sizes
% Following Sanchez-Leon et al.2025 paper: 600 mm^2 ~ 24.5 mm x 24.5 mm
% https://elifesciences.org/reviewed-preprints/100941v2#s2
% 600 mm^2 seems too big, divide by two (for now)
% check in more papers
o.addElectrode(tag = "Lumbar",current = -0.2,type="surface",shape="rectangle", ...
    center=[-1.1819,-9.655,10.129],length=24.5, width = 24.5,thickness=1)

o.run(targetStage=Stage.GETDP,show=false);


