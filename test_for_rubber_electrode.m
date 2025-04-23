% after this script, run rubber_electrode.m
% adjust o.dir and o.ID as needed.

o = EFMouse; % Create a default empty object of the class EFMouse
%o.dir = '/Users/rubensanchez/desktop/EFMouse/4x1Montage_rubber'; % Results and the object (4x1.mat) will be saved here.
%o.ID = '4x1_rubber';  % A name/tag for this simulation.
o.dir = '/Users/rubensanchez/desktop/EFMouse/4x1Montage_copyUp'; % Results and the object (4x1.mat) will be saved here.
o.ID = '4x1_copyUp';  % A name/tag for this simulation.
o.log  = true;  % Create a log file.

o.initialize(overwrite=true);

o.eTag = ["Anterior" "Posterior" "Lateral" "Medial" "Lumbar"];
o.eShape = ["circular" "circular" "circular" "circular" "circular"];
o.eCurrent = [0.05,0.05,0.05,0.05,-0.2];
% for the 4 circular stimulation electrodes
o.eCenter = [-3.56,29.5,5.45;
             -3.43,24.8,6.02;
             -5.5,26.72,3.58; 
              -1.2,27.11,6.29]';
% radius of the electrodes (in mesh space units) (1 mesh unit ~ 1 mm)
o.eRadius = [0.71,0.67,0.64,0.6]';


% for the 1 rectangular return electrode in the lumbar
%o.eXcoord = [-13.4319,11.068]';
%o.eYcoord = [-21.905,2.595]';
%o.eZcoord = [8.5093,9.5093]';

% center coordinate of the craniotomy position (in mesh space)
o.cCenter = [-3.4236,27.1067,5]';
% radius of the creaniotomy (in mesh space units) (1 mesh unit ~ 1 mm)
o.cRadius = 1.5;

% run this before going to rubber_electrode.m
o.run(targetStage=1,show=true)