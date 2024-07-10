classdef EFMouse < handle
    % A class to model electric fields induced by current stimulation in
    % the mouse.
    %
    %
    % Function starting with 'compute' form the pipeline to model the
    % electric field; they are called in the right order by the run function.
    %
    % Function starting with 'plot' visualize results.
    %
    % Function starting with 'analyze' link field measurements with  mouse
    % anatomy.
    %
    % Ruben Sanchez-Romero and Bart Krekelberg
    % Center for Molecular and Behavioral Neuroscience (CMBN)
    % Rutgers Newark
    % March 2024



    properties (SetAccess=public, GetAccess=public)
        eTag (1,:) string      % names for the N stimulation electrodes
        eCurrent (1,:) double  % currents applied to the N electrodes
        eCenter (3,:) double   % 3*N xyz locations of the N electrodes
        eRadius (1,:)         % radii of the N electrodes
        eTissue (1,:) double   % Tissue type of the electrodes
        cCenter (3,:)  double   % XYZ location of the center of the craniotomy
        cRadius  (1,1) double   % Radius of the craniotomy

        dir (1,1) string = tempdir  % Directory where output files will be written
        ID (1,1) string =  "AnID"   % A label for this simulation model

        log (1,1) logical = true; % Create a log file

        % Table of materials with their conductivities. Defaults from ROAST
        conductivityTable dictionary =   dictionary(["gray" "csf" "bone" "skin" "eye" "air"     "conductor" "boundary"],...
            [0.275  1.654  0.01  0.465  0.5    2.5e-14  59e6        59e6]);
        % Map linking a tissue name to an integer ID. Use addTissue to fill.
        tissueLabel  dictionary = dictionary(string([]),[])

        % The mesh is loaded from file (setupBaseModel) then updated in computeMesh and
        % computeBoundaryConditions
        mesh (1,1) struct = struct('node',[],'elem',[],'label',[],'boundary',[],'boundaryLabel',[],'face',[]);
    end

    properties (SetAccess =protected, GetAccess=public)
        stage double {mustBeInteger,mustBeInRange(stage,-1,5)} = -1  % Start at stage = -1
        model

        tissueMaterial  dictionary  = dictionary(string([]),string([]))% Map linking a tissue name to its conductivity. Filled by addTissue
        eArea (1,:) double   % Boundary area for each electrode
        baseModel (1,1) string %Can only be set on construction

    end

    properties (SetAccess =protected, GetAccess=public,Transient)
        % These properties are not saved with the object because they can
        % be read from other saved files.
        field   (:,3) double % Field estimates for each node in o.mesh
        voltage  (:,1) double % Voltage estimate for each node in o.mesh
    end

    properties (Dependent)
        num_electrodes (1,1) double   % Number of stimulation electrodes
        has_craniotomy (1,1) logical  % Does the model have a craniotomy?

    end

    %% Dependent properties
    methods

        function v =get.num_electrodes(o)
            v = numel(o.eTag);
        end

        function v = get.has_craniotomy(o)
            v = ~isempty(o.cCenter);
        end
    end


    %% Construction and Initialization
    methods (Access=public)
        function o= EFMouse(pv)
            % Constructor for an EFMouse object
            %
            % baseModel : the name of the base model to use in a new
            %               simulation ["Alekseichuk"]
            % dir: the folder that contains a previously saved simulation
            % ID: the id of the simulation in the dir folder
            arguments
                pv.baseModel (1,1) string = "Alekseichuk"
                pv.dir (1,1) string = ""
                pv.ID (1,1) string = ""
            end
            if pv.dir ~= ""
                % Load from folder
                o.dir = pv.dir;
                o.ID  = pv.ID;
                load(file(o,"OBJECT"),'o');
                if o.stage >=4
                    % Read the simulation results from file and store in the object
                    o.field = readGetDp(o,type="E");
                    o.voltage = readGetDp(o,type="V");
                end
            else
                % Start clean
                o.baseModel = pv.baseModel;
            end
        end

        function disp(o)
            fprintf('EF Mouse Model (label: %s) in directory %s (stage=%d).\n',o.ID,o.dir,o.stage)
        end

        function initialize(o,pv)
            arguments
                o (1,1) EFMouse
                pv.overwrite (1,1) logical= false
            end
            if o.stage>-1;return;end
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here,'lib/NIfTI_20140122'))
            % make the dir if it does not exist
            if exist(o.dir,'dir')
                if pv.overwrite
                    delete(fullfile(o.dir,'*.*'))
                else
                    error('Folder %s already exists. Load it by specifying dir and id to the constructor, or set overwrite to true to start fresh.\n',o.dir);
                end
            else
                mkdir(o.dir);
            end

            % Setup the base model
            setupBaseModel(o)
            % We are now at stage ==0
            o.field = nan(0,3);
            o.voltage = nan(0,1);
            o.stage  = 0;
        end


        function setupBaseModel(o)
            switch upper(o.baseModel)
                case "ALEKSEICHUK"
                    % load the "clean" mouse mesh (ie. no electrodes or
                    % craniotomy) that Alekseichuk et al created.
                    load('aux_files/EFMouse_mesh_clean.mat','elem','face','node');
                    %Store in the object
                    o.mesh.node = node;
                    o.mesh.elem = elem(1:4,:);
                    o.mesh.label = elem(5,:);
                    o.mesh.face = face;
                    % Map tissue type to labels in the mesh and link them
                    % to entries (with the same name) in the conductivityTable
                    o.addTissue(["gray", "csf"  "bone"  "skin"   "eye"],["gray", "csf"  "bone"  "skin"   "eye"],1:5);
                    % Add tissue types for craniotomy and the skin lesion above
                    % and link them to a material with conductivity.
                    o.addTissue(["craniotomy", "skinremoved"], ["csf" "air"],6:7);
                case "NONE"
                    % Use this if you want to specify the base model "by
                    % hand", without changing any of the class code. This
                    % would require code like the one immediately above in
                    % the user mfile.

                    % case add new base models here
                otherwise
                    error('Unknown baseModel %s (Add it to setupBaseModel?)',o.baseModel)
            end
        end

        function validate(o)
            % Validate that the model specifications meet the requirements
            assert(all(numel(o.eTag) == [numel(o.eCurrent) size(o.eCenter,2) numel(o.eRadius)]),'Each electrode must be assigned a current, a center, and a radius in o.electrodes');
            assert(sum(o.eCurrent)<eps,'Stimulation currents must add up to zero.')
            assert(size(o.cCenter,2)<2,'Only 1 craniotomy can be modeled.');
        end

        function v= file(o,tag)
            % Various files are saved and loaded by different functions, to
            % ensure consistent naming, the filenames are all created here.
            arguments
                o (1,1) EFMouse
                tag (1,1) string {mustBeMember(tag,["OBJECT" "LOG" "E" "V" "TRANS" "DIGIMOUSE" "MESH" "PRO" ...
                    "GETDP" "EMAGNII" "EXNII" "EYNII" "EZNII" ...
                    "ALLEN" "ALLENLABELS"])}
            end
            installDir = fileparts(mfilename("fullpath"));
            switch tag
                case "OBJECT"
                    v = fullfile(o.dir,o.ID + ".mat");
                case "LOG"
                    v = fullfile(o.dir,o.ID + "_logfile.txt");
                case "E"
                    v = fullfile(o.dir,o.ID + "_e.pos");
                case "V"
                    v = fullfile(o.dir,o.ID + "_v.pos");
                case "TRANS"
                    v = fullfile(installDir,"aux_files","transMatrix_ef2Digimouse.mat");
                case "DIGIMOUSE"
                    v = fullfile(installDir,"aux_files","EFMouse_digimouseAtlas.nii.gz");
                case "MESH"
                    v = fullfile(o.dir,o.ID+ ".msh");
                case "PRO"
                    v = fullfile(o.dir,o.ID+ ".pro");
                case "GETDP"
                    % Run getDP
                    str = computer('arch');
                    switch str
                        case 'win64'
                            exe = "getdp.exe";
                        case 'glnxa64'
                            exe = "getdp";
                        case 'maci64'
                            exe = "getdpMac";
                        otherwise
                            error('Unsupported operating system!');
                    end
                    v = fullfile(installDir,"lib","getdp-3.2.0","bin", exe);
                case "EMAGNII"
                    v = fullfile(o.dir,o.ID + "_efm.nii.gz");
                case "EXNII"
                    v = fullfile(o.dir,o.ID + "_efX.nii.gz");
                case "EYNII"
                    v = fullfile(o.dir,o.ID + "_efY.nii.gz");
                case "EZNII"
                    v = fullfile(o.dir,o.ID + "_efZ.nii.gz");
                case "ALLEN"
                    v = fullfile(installDir,"aux_files","EFMouse_allenAtlas.nii.gz");
                case "ALLENLABELS"
                    v = fullfile(installDir,"aux_files","allenAtlas_labels.mat");
                otherwise
                    % cannot happen
            end
        end


        function run(o,pv)
            % RUN runs all stages of the pipeline up to the specified target stage
            % By default the full pipeline is run (all the way to stage 4)
            %   0: Start with a clean mouse mesh
            %   1: Create electrode and craniotomy in the mouse mesh.
            %   2: Compute electrode areas and set boundary conditions .
            %   3: Export mesh and model definitions for getDP
            %   4: Run getDP to solve the electric field model.
            % OPTIONS
            % targetStage : run the pipeline to this stage.
            % clearLog: Clear the log file (if logging is on)
            % startStage: Run the pipeline from this stage.
            % show: plot the mesh and the EF simulation results
            arguments
                o (1,1) EFMouse
                pv.targetStage (1,1) double {mustBeInteger,mustBePositive}  = 4
                pv.clearLog (1,1) logical = false
                pv.startStage (1,1) double = o.stage+1
                pv.show  (1,1) logical = true
            end
            if pv.clearLog
                o.clearLog
            end
            o.stage = min(o.stage,pv.startStage-1);

            if o.stage < pv.targetStage-1
                % Recurse to make sure the pipeline has completed everything up to
                % the stage before the target stage, then run target stage.
                run(o,targetStage = pv.targetStage-1,clearLog = pv.clearLog,show = pv.show);
            elseif pv.targetStage <= o.stage
                % Nothing do
                fprintf('Stage %d already completed.\n',pv.targetStage)
                return;
            end
            %Map stages to compute functions
            tic;
            switch (pv.targetStage)
                case 0
                    % Load the base mesh
                    fprintf('Compute stage %d - %s\n',pv.targetStage,"initialize")
                    o.initialize(overwrite=true);
                case 1
                    % Create the electrodes and craniotomy in the mesh
                    fprintf('Compute stage %d - %s\n',pv.targetStage,"computeMesh")
                    validate(o);
                    computeMesh(o,show=pv.show);
                case 2
                    % Find electrode edges to impose boundary conditions
                    fprintf('Compute stage %d - %s\n',pv.targetStage,"computeBoundary")
                    computeBoundary(o);
                case 3
                    % save the .msh and .pro files
                    fprintf('Compute stage %d - %s\n',pv.targetStage,"export data")
                    saveMesh(o);
                    savePro(o);
                    o.stage = 3;
                case 4
                    % run the GetDP solver to compute electric field
                    % (takes ~30 minutes in a 16G ram mac)
                    fprintf('Compute stage %d - %s\n',pv.targetStage,"computeField")
                    computeField(o,show=pv.show);
                case 5
                    fprintf('Compute stage %d - %s\n',pv.targetStage,"computeVoxelSpace")
                    computeVoxelSpace(o);
                otherwise
                    error('Stage %d??',pv.targetStage)
            end
            save(file(o,"OBJECT"),"o"); % Save so we can pick up later.
            fprintf('Stage %d complete - %.4f seconds \n',o.stage,toc)
        end
    end

    %% Visualization
    methods (Access = public)

        function plotEf(o,pv)
            % PLOTEF visualizes electric field (EF) in mesh space.
            % Focuses on the mouse brain. But can be modified as necessary. See code below.
            %    IN:
            %       percentile: a value from 0 to 100 for the upper threshold of the colormap.
            %           eg. 98 will show EF values < 98th percentile.
            %       orientationQuiver - Set to false to hide the orientation arrows.
            %       tissue:  the name of the tissue to show
            %       type : V (voltage), EMAG (field magnitude), EX, EY,EZ
            %       (field components).
            arguments
                o (1,1) EFMouse
                pv.percentile = 98
                pv.tissue (1,1) string {mustBeMember(pv.tissue,["gray" "csf" "bone"  "skin"   "eye"])} = "gray"
                pv.orientationQuiver (1,1) logical = true
                pv.type (1,1) string {mustBeMember(pv.type,["V" "eMag" "eX" "eY" "eZ"])} = "eMag"
            end

            tic
            fprintf('----Starting plotEf...%s\n',datetime('now'));
            tissueId = o.tissueToLabel(pv.tissue);

            %%% Color limits for the colormap:
            % One of the challenge of the visualization is the distribution
            % of ef values. It has very large outliers usually in areas
            % close to the electrodes. We need to select a percentile as an
            % upper threshold and cut the colormap values there.
            % Otherwise the surface of the brain will look all of one color and
            % we will not be able to visualize the pattern.
            % 98 seems a good trade-off, but you should modify as necessary


            %%
            % Properly, we should only consider in the colormap the values
            % of the tissue we are plotting.
            % This avoid bias in the visualization from other tissue values
            % Get the mesh node indices for the corresponding tissue.
            tiss_node_idx = unique(o.mesh.elem(:,ismember(o.mesh.label,tissueId)));

            % got this colormap from https://jdherman.github.io/colormap/
            load('aux_files/dawn_colormap.mat','dawn_colormap');
            % got this divergent colormap from https://colorbrewer2.org/#type=diverging&scheme=BrBG&n=11
            % and used a Matlab code to expand it to 256 values.
            load('aux_files/brbg_colormap.mat','brbg_colormap');

            switch (pv.type)
                case "eMag"
                    data = sqrt(sum(o.field.^2,2));
                    clims = [0  prctile(data(tiss_node_idx),pv.percentile)];
                    % colormap
                    cm = dawn_colormap;
                case "eX"
                    data = o.field(:,1);
                    clims = prctile(data(tiss_node_idx),[(100-pv.percentile)/2 pv.percentile+(100-pv.percentile)/2]);
                    % divergent colormap to better visualize + and - directions
                    cm = brbg_colormap;
                case "eY"
                    data = o.field(:,2);
                    clims = prctile(data(tiss_node_idx),[(100-pv.percentile)/2 pv.percentile+(100-pv.percentile)/2]);
                    cm = brbg_colormap;
                case "eZ"
                    data = o.field(:,3);
                    clims = prctile(data(tiss_node_idx),[(100-pv.percentile)/2 pv.percentile+(100-pv.percentile)/2]);
                    cm = brbg_colormap;
                otherwise
                    % must be V
                    data = o.voltage;
                    clims = prctile(data(tiss_node_idx),[(100-pv.percentile)/2 pv.percentile+(100-pv.percentile)/2]);
                    cm = dawn_colormap;
            end

            %% Generate the figure
            figure;
            pdeplot3D(o.mesh.node,o.mesh.elem(:,o.mesh.label==tissueId),ColorMapData = data);

            % Note: the following parameters are optional, but we found they allow for
            % a nice visualization of the results.
            % Please modify/comment-out as you consider necessary for your project.

            % define the cutoff (limit) for the colormap values
            clim(clims);
            % define the colorbar
            colormap(cm);
            if ~strcmp(pv.type,'eMag') && ~strcmp(pv.type,'V')
                %cmin = min(data);
                %cmax = max(data);
                % Ensure symmetric range
                cmax = max(abs(clims(1)), abs(clims(2))); 
                cmin = -cmax;
                clim([cmin cmax]);
            end
            % parameters of the colorbar
            cb = colorbar;
            % position and size [left bottom width height]
            cb.Position = [0.82 0.2000 0.0275 0.5000];
            cb.FontSize = 13;
            cb.Box = 'off';
            if strcmp(pv.type,'V')
                % for voltage
                yl = ylabel(cb,'V','Rotation',270);
            else
                % for electric field
                yl = ylabel(cb,'V/m','Rotation',270);
            end
            yl.Position(1) = 4;

            %% Set other parameters of the 3D plot.

            % define title using the montage ID and the tissue
            plot_title = {o.ID + " -- " + pv.tissue,...
                pv.type " < " +  string(num2str(pv.percentile)) + "th percentile"};
            title(plot_title,FontSize=17);
            % rotate for a transversal view (X(left-right)-Y(top-bottom)axes)
            % once plotted, it can be rotated manually with the figure
            % Tools->Rotate 3D option
            view([0,90])
            % fix the position of the figure
            set(gcf,'Position',[440 348 582 449])
            % set the limits for X,Y and Z axis to zoom in on the brain.
            % And also in case the brain is manually rotated (Tools->Rotated 3D)
            % this forces the brain to stay inside this limits
            xlim([-4,6])
            ylim([19,37])
            zlim([0,6])

            % Extra:
            % pdeplot3D has a 3d red orientation marker (plotted as a quiver)
            % helps to navigate if manually rotating the image (Tools->Rotate 3D)
            % Adjust parameters so it can be properly visualized for this brain plot
            % If you do not want it, set pv.orientationQuiver to false
            quiverHandle = findobj(gca, 'type', 'Quiver');
            % Adjust the scale of the quiver plot to make it smaller
            scaleFactor = 0.1;
            set(quiverHandle, 'AutoScaleFactor', scaleFactor);
            % x,y,zlim values can be used as reference to set the new values
            newYData = quiverHandle.YData + 90;
            newZData = quiverHandle.ZData + 40;
            newXData = quiverHandle.XData + 47;
            % Update the data of the quiver to move it to the new position
            set(quiverHandle, 'YData', newYData,'ZData',newZData,'XData',newXData);
            % modify line width
            set(quiverHandle,'LineWidth', 1);
            % modify arrowhead size
            set(quiverHandle,'MaxHeadSize',0.3);
            % add labels X,Y,Z to the red quiver
            text(6, 20, 4, 'X'); % try "right"
            text(3, 23.1, 4, 'Y'); %try "anterior"
            text(3, 20, 6.9, 'Z'); % try "superior"

            if ~pv.orientationQuiver
                delete(findobj(gca,'type','Quiver'));
            end

            fprintf('----plotEf elapsed time: %.4f seconds\n', toc);
        end

        function plotMesh(o,pv)
            % PLOTMESH visualizes the montage in mesh space.
            % It zooms-in on the brain even though return electrodes could be
            % in the lumbar or some other posterior or inferior parts of the mouse.
            % To visualize other sections of the body the user can modify the
            % code or use the graphical display tools to zoom-out and navigate.
            % Color code:
            %   Craniotomy:
            %       skin removal plotted: Green
            %       skull removal plotted: Yellow
            %   Electrodes with + current are plotted in Red
            %   Electrodes with - current are plotted in Blue
            %   Roi plotted: Magenta
            %
            arguments
                o (1,1) EFMouse
                pv.xlim (1,2) double = [-4,6]
                pv.ylim (1,2) double =  [19,37] %[23,45] %[19,37]
                pv.zlim (1,2) double = [0,6]
                pv.roi (1,:) double = []
            end

            tic
            clf;
            % Plot the soft tissue of the mesh
            [tf,ix] = meshHasTissue(o,"gray");
            if tf
                % ix=ix(1:6:end);
                pdeplot3D(o.mesh.node,o.mesh.elem(:,ix),'FaceColor','white','FaceAlpha',0.001);
                hold on;
            end

            % plot the craniotomy skin removal (if any)
            [tf,ix] = meshHasTissue(o,"skinremoved");
            if tf
                pdeplot3D(o.mesh.node,o.mesh.elem(:,ix),'FaceColor','green','EdgeColor','green');
                hold on;
            end

            % plot the craniotomy skull removal (if any)
            [tf,ix] = meshHasTissue(o,"craniotomy");
            if tf
                pdeplot3D(o.mesh.node,o.mesh.elem(:,ix),'FaceColor','yellow','EdgeColor','yellow');
                hold on;
            end

            % plot the electrodes (+ current: red) and (- current: blue)
            for i = 1:o.num_electrodes
                [tf,ix] = meshHasTissue(o,o.eTag(i));
                if tf
                    current = o.eCurrent(i);
                    if current > 0 % positive current
                        thisColor ='red';
                    elseif current < 0 % negative current
                        thisColor ='blue';
                    end
                    pdeplot3D(o.mesh.node,o.mesh.elem(:,ix),'FaceColor',thisColor,'EdgeColor',thisColor);
                    hold on;
                end
            end
            
            % plot the roi, when analyzing roi-level electric field
            if ~isempty(pv.roi)
                pdeplot3D(o.mesh.node,o.mesh.elem(:,pv.roi),'FaceColor','magenta','EdgeColor','magenta');
                hold on;
            end

            %% Parameters of the plot

            % define title using the montage ID
            %plot_title = o.ID;
            %title(plot_title,FontSize=17);
            % change the position of the title, it overlaps with the mouse head
            %title_ax = get(gca,'Title');
            %title_ax.Position = [8 45.1127 3.0000];

            %if ~isempty(o.eTag)
            %    % Beta: still thinking which is the best to display this info
            %    % in a table? next to the electrodes?
            %    % add the tags of the electrodes and the currents (code from chatGPT)
            %    % TODO: confirm the units used in ROAST
            %    list_electrode_tags = sprintf('%s\n', o.eTag(:));
            %    text(12, 35, 3, list_electrode_tags, 'FontSize', 13, 'FontWeight', 'bold');
            %    list_electrode_currents = sprintf('%.2f mA\n',o.eCurrent);
            %    text(8, 35, 3, list_electrode_currents, 'FontSize', 13, 'FontWeight', 'bold');
            %end
            % rotate for a transversal view (X(left-right)-Y(anterior-posterior)axes)
            % once plotted, the body can be rotated manually with Tools->Rotate 3D option
            view([0,90])
            % fix the position of the figure
            set(gcf,'Position',[440 348 582 449])
            % Set the limits for X,Y and Z axis to zoom in on the head.
            % This forces the brain to stay inside this limits. Change as needed.
            xlim(pv.xlim)
            ylim(pv.ylim)
            zlim(pv.zlim)

        end
    end

    %% Functions to add or search for tissue types
    methods (Access= public)
        function v = labelToTissue(o,id)
            % Given a tissue type id (a number), return its tissue label (
            % a string)
            arguments
                o (1,1) EFMouse
                id (1,:) double {mustBeInteger}
            end
            allIds = o.tissueLabel.values;
            [tf,loc] = ismember(id,allIds);
            v = repmat("Tissue not defined",size(id));
            if any(tf)
                allKeys= o.tissueLabel.keys;
                v(tf) = string(allKeys(loc));
            end
        end
        function v = tissueToLabel(o,name)
            % Given a tissue name (a string) return its tissue ID (
            % a number)
            arguments
                o (1,1) EFMouse
                name (1,:) string
            end
            isAKey = isKey(o.tissueLabel,name);
            v =NaN(size(name));
            for i= find(isAKey)
                v(i) = o.tissueLabel(name(i));
            end

        end
        function  v = addTissue(o,name,material,label)
            % Add a tissue type to the model and return the ID assigned to
            % the type.
            % name : The name of the tissue
            % material:  the kind of material the tissue is made of. This
            % should be one of the elements in the conductivityTable.
            % value: A unique label (integer) that identifies this
            % tissue in the mesh.
            arguments
                o (1,1) EFMouse
                name (1,:) string
                material (1,:) string % Material determines the conductivity
                label (1,:) double {mustBePositive,mustBeInteger}=[]
            end
            if isempty(label)
                label = o.tissueLabel.numEntries +(1:numel(name));
            end
            assert(numel(name)==numel(label),"The number of names should match the number of values ")

            alreadyDefined = isKey(o.tissueLabel,name);
            % Add to the dictionary
            o.tissueLabel(name(~alreadyDefined)) = label(~alreadyDefined);
            v= tissueToLabel(o,name);
            if numel(material)==1
                material =repmat(material,[1 numel(name)]);
            end
            o.tissueMaterial(name) = material;
        end


        function [tf,ix] = meshHasTissue(o,nameOrId)
            if isstring(nameOrId)
                nameOrId = o.tissueToLabel(nameOrId);
            end
            ix = find(o.mesh.label == nameOrId);
            tf = ~isempty(ix);
        end
    end

    %% Helper functions
    methods (Access =protected)
        function startLog(o)
            if o.log
                diary(file(o,"LOG"));
            end
        end

        function stopLog(o)
            if o.log
                diary('off')
            end
        end
        function clearLog(o)
            if o.log && exist(file(o,"LOG"),"file")
                del(file(o,"LOG"))
            end
        end
        function [data,nodeNr] = readGetDp(o,pv)
            % Read .pos files that contain the output of GetDP. The first
            % number is the number of nodes, followed by nrNodes lines
            % representing the node number (column 1) and then the data
            % (e.g. 1 column for voltage , 3 columns for field).
            %
            % type: "V" for voltage, "E" for electric field.
            arguments
                o (1,1) EFMouse
                pv.type (1,1) string {mustBeMember(pv.type,["V" "E"])}
            end

            fname = file(o,pv.type);
            fid = fopen(fname);
            tmp = fscanf(fid,'%f');
            fclose(fid);
            nrNodes= tmp(1);
            data = reshape(tmp(2:end),[],nrNodes)';
            nodeNr = data(:,1);
            data(:,1)=[];
        end

    end

    %% Core pipeline
    methods (Access=protected)



        function saveMesh(o)
            % SAVEMESH saves the mesh in *.msh format. This is the
            % file used by the electric field modeling GetDP solver.

            % Uses code from ROAST

            % start logging
            o.startLog
            tic
            fprintf('----Starting saveMesh...%s\n',datetime('now'));
            % save the mesh using iso2mesh function
            meshFile  =file(o,"MESH");
            fid = fopen(meshFile,'wt');
            fprintf (fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');
            % Write the nodes
            nrNodes= size(o.mesh.node,2);
            fprintf (fid, '$Nodes\n');
            fprintf (fid, '%d\n', nrNodes);
            buffer = [ (1:nrNodes)' o.mesh.node']';
            fprintf (fid, '%d %10.10f %10.10f %10.10f\n', buffer);
            fprintf (fid, '$EndNodes\n');

            % Write the tetrahedra [4 2] means "tetrahedron" and "2
            % tags"  (in our case both tags are the same; the label)
            nrTetrahedra = size(o.mesh.elem,2);
            nrBoundary  = size(o.mesh.boundary,2);
            fprintf (fid, '$Elements\n');
            fprintf (fid, '%d\n',nrTetrahedra+nrBoundary);
            buffer = [(1:nrTetrahedra)' repmat([4 2],nrTetrahedra,1) repmat(o.mesh.label',[1 2]) o.mesh.elem']';
            fmt = [repmat('%d ',[1 9] ) '\n'];
            fprintf (fid, fmt, buffer);

            % Write the boundary triangles - [2 2] means "triangle" and "2
            % tags"  (in our case both tags are the same; the label)
            buffer = [nrTetrahedra+(1:nrBoundary)' repmat([2 2],nrBoundary,1) repmat(o.mesh.boundaryLabel',[1 2]) o.mesh.boundary']';
            fmt = [repmat('%d ',[1 8] ) '\n'];
            fprintf (fid, fmt, buffer);
            fprintf (fid, '$EndElements\n');
            fclose(fid);
            fprintf('----saveMesh elapsed time: %.4f seconds\n\n', toc);
            o.stopLog
        end

        function proFile = savePro(o)
            % PROFILE creates a text file (*.pro) with the running call for getDP
            % (This defines which pde model need to be computed).

            % start logging
            o.startLog
            tic
            fprintf('----Starting savePro...%s\n',datetime('now'));
            
            % Define file names for function output
            % a text file that will contain input parameters for the FEM solver
            proFile = file(o,"PRO");
            % file to save the voltage results of the model
            % do not need to include the path, only the name
            % the GetDP solver finds the path from the .pro file
            output_v = file(o,"V");
            % file to save the electric field results of the model
            % same as above, no need to include the path
            output_e = file(o,"E");

            %%
            % assign these variables names to keep ROAST convention
            % extract the name of the electrodes
            % TODO use electrode names  elecName = o.eTags';

            fid = fopen(proFile,'w');
            fprintf(fid,'/* \n .pro file created by EFMouse on %s \n ID: %s \n Dir: %s \n */\n\n',datetime("now"),o.ID,o.dir);
            %% Define the tissues (Region) in GetDP format
            fprintf(fid,'Group {\n\n');
            tissues= o.tissueLabel.keys;
            tissuesWithoutBoundaries= tissues(~startsWith(tissues,'boundary'));
            nrTissWithoutBoundaries= numel(tissuesWithoutBoundaries);
            nrTiss = numel(tissues);
            for k= 1:nrTiss
                fprintf(fid,'%s = Region[%d];\n', tissues{k},o.tissueToLabel(tissues{k}));
            end
            fprintf(fid,'DomainC = Region[{');
            fprintf(fid,'%s',strjoin(tissuesWithoutBoundaries,','));
            fprintf(fid,'}];\n');

            fprintf(fid,'AllDomain = Region[{');
            fprintf(fid,'%s',strjoin(tissues,','));
            fprintf(fid,'}];\n}\n\n');
            %% Define conductivities for each of the regions
            fprintf(fid,'Function {\n\n');
            for k= 1:nrTissWithoutBoundaries
                fprintf(fid,'sigma[%s] = %g;\n',tissuesWithoutBoundaries{k},o.conductivityTable(o.tissueMaterial(tissuesWithoutBoundaries{k})));
            end
            %% Define the currents for the electrode surfaces (boundary elements)
            for i=1:o.num_electrodes
                fprintf(fid,'du_dn%d[] = %f;\n',i,1000*o.eCurrent(i)/o.eArea(i));
            end
            fprintf(fid,'}\n\n');

            %% Define the electrostatic model parameters
            fprintf(fid,'Jacobian {\n');
            fprintf(fid,'  { Name Vol ;\n');
            fprintf(fid,'    Case {\n');
            fprintf(fid,'      { Region All ; Jacobian Vol ; }\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'  { Name Sur ;\n');
            fprintf(fid,'    Case {\n');
            fprintf(fid,'      { Region All ; Jacobian Sur ; }\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'}\n\n\n');

            fprintf(fid,'Integration {\n');
            fprintf(fid,'  { Name GradGrad ;\n');
            fprintf(fid,'    Case { {Type Gauss ;\n');
            fprintf(fid,'            Case { { GeoElement Triangle    ; NumberOfPoints  3 ; }\n');
            fprintf(fid,'                   { GeoElement Quadrangle  ; NumberOfPoints  4 ; }\n');
            fprintf(fid,'                   { GeoElement Tetrahedron ; NumberOfPoints  4 ; }\n');
            fprintf(fid,'                   { GeoElement Hexahedron  ; NumberOfPoints  6 ; }\n');
            fprintf(fid,'                   { GeoElement Prism       ; NumberOfPoints  9 ; } }\n');
            fprintf(fid,'           }\n');
            fprintf(fid,'         }\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'}\n\n\n');

            fprintf(fid,'FunctionSpace {\n');
            fprintf(fid,'  { Name Hgrad_v_Ele; Type Form0;\n');
            fprintf(fid,'    BasisFunction {\n');
            fprintf(fid,'      // v = v  s   ,  for all nodes\n');
            fprintf(fid,'      //      n  n\n');
            fprintf(fid,'      { Name sn; NameOfCoef vn; Function BF_Node;\n');
            fprintf(fid,'        Support AllDomain; Entity NodesOf[ All ]; }\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'}\n\n\n');

            fprintf(fid,'Formulation {\n');
            fprintf(fid,'  { Name Electrostatics_v; Type FemEquation;\n');
            fprintf(fid,'    Quantity {\n');
            fprintf(fid,'      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'    Equation {\n');
            fprintf(fid,'      Galerkin { [ sigma[] * Dof{d v} , {d v} ]; In DomainC; \n');
            fprintf(fid,'      Jacobian Vol; Integration GradGrad; }\n');
            % Currents.
            for i=1:o.num_electrodes
                fprintf(fid,'      Galerkin{ [ -du_dn%d[], {v} ]; In %s ;',i,"boundary"+ o.eTag(i));
                fprintf(fid,'                 Jacobian Sur; Integration GradGrad;}\n');
            end

            fprintf(fid,'    }\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'}\n\n\n');

            fprintf(fid,'Resolution {\n');
            fprintf(fid,'  { Name EleSta_v;\n');
            fprintf(fid,'    System {\n');
            fprintf(fid,'      { Name Sys_Ele; NameOfFormulation Electrostatics_v; }\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'    Operation { \n');
            fprintf(fid,'      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'}\n\n\n');

            fprintf(fid,'PostProcessing {\n');
            fprintf(fid,'  { Name EleSta_v; NameOfFormulation Electrostatics_v;\n');
            fprintf(fid,'    Quantity {\n');
            fprintf(fid,'      { Name v; \n');
            fprintf(fid,'        Value { \n');
            fprintf(fid,'          Local { [ {v} ]; In AllDomain; Jacobian Vol; } \n');
            fprintf(fid,'        }\n');
            fprintf(fid,'      }\n');
            fprintf(fid,'      { Name e; \n');
            fprintf(fid,'        Value { \n');
            fprintf(fid,'          Local { [ -{d v} ]; In AllDomain; Jacobian Vol; }\n');
            fprintf(fid,'        }\n');
            fprintf(fid,'      }\n');
            fprintf(fid,'    }\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'}\n');

            fprintf(fid,'PostOperation {\n');
            fprintf(fid,'{ Name Map; NameOfPostProcessing EleSta_v;\n');
            fprintf(fid,'   Operation {\n');
            fprintf(fid,'     Print [ v, OnElementsOf DomainC, File "%s", Format NodeTable ];\n',output_v);
            fprintf(fid,'     Print [ e, OnElementsOf DomainC, Smoothing, File "%s", Format NodeTable ];\n',output_e);
            fprintf(fid,'   }\n');
            fprintf(fid,'}\n\n\n');
            fprintf(fid,'}\n');
            % close the file
            fclose(fid);

            fprintf('----savePro elapsed time: %.4f seconds\n\n', toc);
            o.stopLog
        end


        function computeMesh(o,pv)
            % COMPUTEMESH creates a user-defined stimulation electrodes and 
            % craniotomy on mesh space.
            %
            % show: show the resulting mesh
            arguments
                o (1,1) EFMouse
                pv.show = false;
            end

            o.startLog;
            disp(o)

            % create a pde model using the mesh, this allows us to use the
            % Matlab functions that create the electrodes and the craniotomy
            o.model = createpde();
            geometryFromMesh(o.model,o.mesh.node,o.mesh.elem,o.mesh.label);


            %% If defined, create the craniotomy: we are just modeling one craniotomy
            % craniotomy may or may not be defined, depending on the experiment
            % simulated.
            % Check if craniotomy information is available
            if o.has_craniotomy
                % do craniotomy and skin removal if information is available
                craniotomy = findElements(o.model.Mesh,'radius',...
                    o.cCenter,...
                    o.cRadius);
                % Check the craniotomy includes bone(skull) and skin.
                % We do not care if it touches grey matter (1) or CSF (2),
                % since we will not modify those labels.
                tissue_touched = unique(o.mesh.label(craniotomy));
                fprintf('-Creating craniotomy-\n')
                fprintf(' Craniotomy: touching tissue:\n');
                for t = 1:numel(tissue_touched)
                    tiss = tissue_touched(t);
                    elem_tiss_touched = sum(o.mesh.label(craniotomy) == tiss);
                    fprintf('   %s: num elements = %d\n', o.labelToTissue(tiss),elem_tiss_touched);
                end
                % Skin removal modeling:
                % Change skin label on craniotomy to label (skinremoved)
                thisId = addTissue(o,"skinremoved","air");
                change = craniotomy(o.mesh.label(craniotomy)==o.tissueToLabel("skin"));
                o.mesh.label(change) = thisId;
                % Skull removal modeling:
                % Change skull/bone label on craniotomy to label
                % then assign conductivity of CSF.
                thisId = addTissue(o,"craniotomy","csf");
                change = craniotomy(o.mesh.label(craniotomy)==o.tissueToLabel('bone'));
                o.mesh.label(change) = thisId;
                fprintf(' Craniotomy only removes skin and bone.\n')
            end

            %% define the n electrodes
            o.eTissue=nan(1,o.num_electrodes); % Tisse type that was used to create the electrodes.
            fprintf('-Creating %d electrodes-\n',o.num_electrodes);
            for i = 1:o.num_electrodes
                % build the electrode by creating a sphere with user-input center
                % and radius.
                % https://www.mathworks.com/help/pde/ug/pde.femesh.findelements.html
                % ("electrode" contains the mesh elements comprising the electrode.)
                electrode = findElements(o.model.Mesh,'radius',...
                    o.eCenter(:,i),...
                    o.eRadius(i));

                % the assigned electrode label in the mesh
                thisId = o.addTissue(o.eTag(i), "conductor");% Add electrode labels at the end of the "tissue" list.

                % Type of tissue "touched" by the electrode:
                % this is a way to assess if our electrode is in the tissue we want
                % if not, we need to modify center and radius.
                tissue_touched = unique(o.mesh.label(electrode));
                elem_tiss_touched  = zeros(1,numel(tissue_touched));

                fprintf(' Electrode: %s: touching tissue:\n',o.eTag(i));
                for t = 1:numel(tissue_touched)
                    tiss = tissue_touched(t);
                    elem_tiss_touched(1,t) = sum(o.mesh.label(electrode) == tiss);
                    fprintf('   %s: num elements = %d\n', o.labelToTissue(tiss),elem_tiss_touched(1,t));
                    % change the tissue labels of the elements to the new electrode tag
                end
                % create electrode only in the max touched tissue
                [~,idx] = max(elem_tiss_touched);
                tiss_elec = tissue_touched(idx);
                fprintf('   Creating electrode %s in max touched tissue: %s.\n',o.eTag(i),o.labelToTissue(tiss_elec))
                change = electrode(o.mesh.label(electrode)==tiss_elec);
                o.mesh.label(change) = thisId;
                % Save the tissue where the electrode was inserted, this will be used
                % by computeBoundary to define the boundary conditions.
                o.eTissue(i) = tiss_elec;
            end

            o.stage = 1;

            if pv.show
                % visualize the resulting mesh
                % Helpful to confirm that electrodes and craniotomy are
                % placed correctly. If not, modify the object and rerun.
                plotMesh(o);
                drawnow;
            end

            % stop logging
            o.stopLog;
          
        end


        function computeBoundary(o)
            % COMPUTEBOUNDARY computes electrode areas and sets boundary conditions for
            % the mesh. Based on the ROAST function prepareForGetDP.m

            o.startLog
            % preallocate
            o.eArea= zeros(o.num_electrodes,1);
            %%
            % turn off the warning
            % for when points are provided outside the bounds of a triangulation.
            warnStt =warning('query');
            warning('off','MATLAB:TriRep:PtsNotInTriWarnId');
            %Clean the boundary element in the mesh
            o.mesh.boundary = [];
            o.mesh.boundaryLabel = [];
            % loop through each electrode
            for i = 1:o.num_electrodes
                % Find electrode nodes and nodes with the same tissue type
                [~,indNode_elecElm] = o.meshHasTissue(o.eTag(i));
                [~,indNode_tissElm] = o.meshHasTissue(o.eTissue(i)); % get the tissue where the electrode was inserted in computeMesh
                % Error if the electrode does not have any nodes
                assert(~isempty(indNode_elecElm),sprintf('Electrode %d was not meshed properly. Reasons may be: 1) electrode size is too small so the mesher cannot capture it; 2) mesh resolution is not high enough. Consider using bigger electrodes or increasing the mesh resolution by specifying the mesh options.',i));

                % In our mouse model no gel is used.
                % instead of gel as in ROAST, use the tissue where the electrode was
                % inserted to setup the boundary conditions

                [~,verts_tiss] = freeBoundary(TriRep(o.mesh.elem(:,indNode_tissElm)',o.mesh.node'));
                [faces_elec,verts_elec] = freeBoundary(TriRep(o.mesh.elem(:,indNode_elecElm)',o.mesh.node'));

                % Find the faces between the electrode and the surrounding tissue
                [~,iE,~] = intersect(verts_elec,verts_tiss,'rows');
                tempTag = ismember(faces_elec,iE);
                faces_elecOuter = faces_elec(~(sum(tempTag,2)==3),:);
                [~,Loc] = ismember(verts_elec,o.mesh.node','rows');
                boundaryElements = Loc(faces_elecOuter);
                % calculate the surface area
                a = (verts_elec(faces_elecOuter(:, 2),:) - verts_elec(faces_elecOuter(:, 1),:));
                b = (verts_elec(faces_elecOuter(:, 3),:) - verts_elec(faces_elecOuter(:, 1),:));
                c = cross(a, b, 2);
                o.eArea(i) = sum(0.5*sqrt(sum(c.^2, 2)));
                assert(o.eArea(i) > 0,'Error: electrode %d area needs to be > 0. Make the radius of the electrode larger and try again.',i);

                thisLabel  = o.addTissue("boundary" + o.eTag(i),"boundary");
                o.mesh.boundary = [o.mesh.boundary boundaryElements'];
                o.mesh.boundaryLabel = [o.mesh.boundaryLabel repmat(thisLabel,[1  size(boundaryElements,1)])];
            end
            o.stage = 2; % This completes stage 2.
            
            % Print the electrode boundary sizes
            T= table(o.eArea', 'RowNames',cellstr(o.eTag)','VariableNames',{'Electrode Area'});
            disp(T)
            o.stopLog
            warning(warnStt);
        end

        function computeField(o,pv)
            % solve solves the electric field model using GetDP. This is an
            % adaptation of ROAST solveByGetDP.m function.
            % See: https://getdp.info/ for info about the FEM solver.
            %   OUT: (*_.pro) txt input file to run GetDP (run here).
            %        (*_e.pos) txt file with X,Y,Z components for each mesh node
            %               electric field.
            %        (*_v.pos) txt file with voltage for each mesh node.
            arguments
                o (1,1) EFMouse
                pv.show = false;
            end
            o.startLog
            % Sanity check
            assert(all(isKey(o.tissueMaterial,o.tissueLabel.keys)),"All tissues in tissueLabel must have an entry in tissueType to link it to the impedanceTable");
            assert(o.tissueLabel.numEntries==numel(unique(o.tissueLabel.values)),"All tissues in tissueLabel must have a unique label number")
            proFile = file(o,"PRO");
            assert(exist(proFile,"file"),"Could not find %s",proFile);

            exe  = file(o,"GETDP");
            cmd = sprintf('"%s" "%s" -solve EleSta_v -msh "%s" -pos Map',exe ,proFile, fullfile(o.dir,o.ID+ ".msh"));
            try
                [status,cmdout] = system(cmd,'-echo');
            catch
                status=1;cmdout = "Call to getdp generated an error";
            end
            assert(status==0,'getDP failed (%s). Please check error messages above.',cmdout);
            % Read the results and store in the object
            o.field = readGetDp(o,type="E");
            o.voltage = readGetDp(o,type="V");
            o.stopLog
            o.stage = 4;
            if pv.show
                plotEf(o)
            end
        end

        function computeVoxelSpace(o,pv)
            % COMPUTEVOXELSPACE transforms EF model results from mesh space to
            % volume-voxel (3D matrix) space.
            % Follows some code from ROAST function postGetDP.m.
            % The voxel-space is aligned to Digimouse-Allen stereotaxic-Paxinos.
            % (Important!: Requires FSL for the alignment.)
            %  tissue: the name of the tissue type for which the field values will be exported. (Defaults to gray)
            %   OUT: (*.nii.gz) 4 nifti files for EF: x,y,z components and magnitude.
            arguments
                o (1,1) EFMouse
                pv.tissue (1,:) string = "gray"
                pv.keepUnaligned (1,1) logical = true
            end
            tic
            o.startLog
            fprintf('----Starting computeVoxelSpace...%s\n',datetime('now'));

            % extract the node indices corresponding to the Grey matter elements
            tissue_node_idx= unique(o.mesh.elem(:,ismember(o.mesh.label,o.tissueToLabel(pv.tissue))));
            nodeXYZ = o.mesh.node(:,tissue_node_idx);

            %% mesh to voxel space transformation
            % Convert mesh coordinates to voxel coordinates for
            % interpolation into voxel space.

            % Pixdim from aux_files/EFMouse_digimouseAtlas.nii.gz,
            % Information is in the header of that nifti file
            % It is in milimeters (0.1 mm)
            % TODO: opent the nifti and extract the info. Instead of hard coding!
            % NOTE: Not sure this step is necessary
            vox_dim = 0.1;
            nodeXYZ = nodeXYZ./vox_dim;

            % Get min and max coordinates for the nodes.
            % This is necessary to create a grid that covers the full mesh grey matter
            minXYZ = min(nodeXYZ,[],2);
            maxXYZ = max(nodeXYZ,[],2);


            % Set the "Grid Resolution" so that it matches the number of voxels in
            % the aux_files/EFMouse_digimouseAtlas.nii.gz
            % This information can be extracted from the header of that nifti file
            % TODO: opent the nifti and extract the info. Instead of hard coding!
            multiple = 1; % change this if want to increase the resolution
            resolXYZ = (maxXYZ-minXYZ)./(multiple*[108 191 86]'-1);

            % create the grid with the appropriate size (full coverage of the mesh
            % brain) and the appropriate resolution (number of voxels in the digimouseAtlas)
            % This grid contains 3D indices for each voxel in the voxel-space
            [xi,yi,zi] = ndgrid(minXYZ(1):resolXYZ(1):maxXYZ(1),minXYZ(2):resolXYZ(2):maxXYZ(2),minXYZ(3):resolXYZ(3):maxXYZ(3));

            %% Interpolate EF mesh values to the grid points of voxel-space
            % this are the computations that take more time in the function

            % Array where the interpolated voxel-space values will be saved
            % the last dimension will save EF values for X,Y and Z direction
            ef = zeros([size(xi) 3]);

            % Computes a linear interpolant function F
            % based on the EF X direction values of the grey matter nodes
            for dim=1:3
                F = TriScatteredInterp(nodeXYZ', o.field(tissue_node_idx,dim));
                % Apply the interpolant F for the xi,yi,zi voxel-space coordinates
                ef(:,:,:,dim) = F(xi,yi,zi);
            end

            fprintf('----Exporting to NIFTI volumes \n');
            %% create nifti files and align to digimouse-atlas
            % Use auxiliary function to create nifti files with the proper header
            % and recentering (stereotaxic-Paxinos layout approximation)
            for dim=0:3
                % Extract the relevant component
                if dim==0
                    % Compute the EF magnitude as sqrt(x^2 + y^2 + z^2), using the
                    % voxel-space interpolated values of ef
                    thisE = sqrt(sum(ef.^2,4));
                    fname = file(o,"EMAGNII");

                else
                    thisE = ef(:,:,:,dim); 
                    fname = file(o,"E" + char('X'+dim-1) + "NII");
                end
                fname = strrep(fname,'.nii','_unaligned.nii');
                % set NaNs to zero (these come from the interpolant function F)
                thisE(isnan(thisE))=0;
                % Create and save nifti
                niiThisE= EFMouse.niftiCreator(thisE);
                save_nii(niiThisE ,char(fname));
                % Use FSL to align the nifti to  digimouseAtlas.
                align2DigimouseAtlas_fsl(o,char(fname));
                if ~pv.keepUnaligned
                    delete(fname);  % Delete unaligned file
                end
            end
            fprintf('----computeVoxelSpace elapsed time: %.4f seconds\n', toc);
            o.stopLog
            o.stage =5;
        end
    end




    %% Analysis functions
    methods (Access = public)

        function analyzeTissue(o)
            % ANALYZETISSUE reports electric field (EF) x,y,z components and
            % magnitude summary statistics (mean, median, std.dev, min, max)
            % for all tissues in the Digimouse.
            arguments
                o (1,1) EFMouse
            end
            tic
            fprintf('----Starting analyzeTissue...%s\n',datetime('now'));
            %% Read the results
            ef = o.readGetDp(type ="E");
            keyDic = keys(o.tissueLabel);
            n_tissue = 5;
            if o.has_craniotomy
                n_tissue = n_tissue + 2;
            end
            for i = 1:n_tissue
                % Select only the corresponding tissue type
                tissue = keyDic(i);
                tissue_node_idx= unique(o.mesh.elem(:,ismember(o.mesh.label,o.tissueToLabel(tissue))));
                ef_tissue = ef(tissue_node_idx,:);
                S = EFMouse.summary(ef_tissue);
                fprintf('Electric field summary statistics for %s tissue\n',tissue);
                disp(S)
                % Compute homogeneity
                H = EFMouse.homogeneity(ef_tissue);
                fprintf('[Homogeneity ranges from 0 to 1]\n');
                fprintf('Homogeneity = %.4f, for ef_norm_mean: (%.3f %.3f %.3f)\n\n',H(1),H(2),H(3),H(4));

            end

        end


        function [S,F,H] = analyzeRoi(o,roi,pv)
            % ANALYZEROI reports electric field (EF) x,y,z components and
            % magnitude summary statistics (mean, median, std.dev, min, max) for a
            % user defined region of interest in mesh space.
            % It also computes a focality measure relative to the given roi
            %   IN: o: the model
            %       roi: a structure with information about the roi. See below. Needs
            %       to be created before running this function. A spherical and a
            %       box-shaped option.
            %   OUT: a summary table (Matlab table object).
            %        a mesh figure showing the roi in the mouse brain.
            %        (including electrodes and craniotomy as reference.)
            %
            % eg1.  roi.shape = 'spherical';
            %       roi.center = [-1.8276,28.7571,3.30];
            %       roi.radius = 2;
            %       summary_table = analyzeRoi(o,roi);
            %
            % eg2.  roi.shape = 'box';
            %       roi.dim = [[-2.77 -1.29];[27.55 29.64];[3.68 4.68]]'
            %       summary_table = analyzeRoi(o,roi);
            arguments
                o (1,1) EFMouse
                roi (1,1) struct
                pv.tissue (1,:) string = "gray"
                pv.plot (1,1) logical = true
                pv.foc_threshold = 75;
                pv.foc_percentile_max = 99.9;
            end

            tic
            fprintf('----Starting analyzeRoi...%s\n',datetime('now'));

            if strcmp(roi.shape,'spherical')
                roi_idx = findElements(o.model.Mesh,'radius',roi.center,roi.radius);
            elseif strcmp(roi.shape,'box')
                roi_idx = findElements(o.model.Mesh,'box',roi.dim(:,1),roi.dim(:,2),roi.dim(:,3));
            end
            fprintf('Electric field summary statistics for a %s roi in %s\n',roi.shape,pv.tissue);

            %% Read the results
            ef = o.readGetDp(type ="E");
            % Select only the requested tissue type
            roi_idx = roi_idx(ismember(o.mesh.label(roi_idx),o.tissueToLabel(pv.tissue)));
            % these roi_idx are tetrahedra elements, need to get 
            % the corresponding nodes
            roi_node_idx = unique(o.mesh.elem(:,roi_idx));
            % extract the corresponding electric field information
            roi_ef= ef(roi_node_idx,:);
            S = EFMouse.summary(roi_ef);
            disp(S)
            if pv.plot
                plotMesh(o,roi=roi_idx);
            end
            % Compute the relative focality metric
            % need to get the complement of the roi (only for the gray matter)
            % this is at the node level.
            gray_node_idx = unique(o.mesh.elem(:,ismember(o.mesh.label,o.tissueToLabel('gray'))));
            % get the elements in gray that are not in roi (complement)
            complement_node_idx = setdiff(gray_node_idx,roi_node_idx);
            complement_ef= ef(complement_node_idx,:);
            F = EFMouse.focality(roi_ef,complement_ef,threshold = pv.foc_threshold, percentile_max = pv.foc_percentile_max);
            fprintf('Relative focality = %.4f%% of %d reference nodes (eMag > %.2f%% of the ROI max (%.2fth percentile))\n',F(1)*100,F(2),pv.foc_threshold,pv.foc_percentile_max);
            % Compute homogeneity
            H = EFMouse.homogeneity(roi_ef);
            fprintf('[Homogeneity ranges from 0 to 1]\n');
            fprintf('Homogeneity = %.4f, for ef_norm_mean: (%.3f %.3f %.3f)\n',H(1),H(2),H(3),H(4));
        end

        function [S,F,H] = analyzeAtlas(o,area,pv)
            % analyzeAtlas reports electric field (EF) x,y,z components and
            % magnitude summary statistics (mean, median, std.dev, min, max) for a
            % chosen hemisphere and Allen atlas area of interest.
            % Requires the NIfTI files output by computeVoxelSpace
            % (stage 5 of the compute pipeline)
            %   IN: o - the EFMouse object
            %       hemisphere: a string "left", "right" or "both".
            %       area: a vector of strings representing Allen atlas areas
            %       (see https://scalablebrainatlas.incf.org/composer/?template=ABA_v3).
            %   OUT: a summary table (Matlab table object)
            %
            % eg. T= analyzeAtlas(o,"Visual areas",hemisphere="LEFT");
            arguments
                o (1,1) EFMouse
                area (1,:) string
                pv.hemisphere (1,1) string {mustBeMember(pv.hemisphere,["left","right","both"])} = "both"
                pv.foc_reference (1,1) string = "Isocortex"
                pv.foc_threshold = 75;
                pv.foc_percentile_max = 99.9;

            end
            tic
            fprintf('----Starting analyzeAtlas...%s\n',datetime('now'));

            %% Find areas and Atlas ids that belong to the requested area
            [allenNames,allenIds] = queryAllen(o,area);
            if allenNames=="";S=table;return;end

            %% Load the voxel matrices from the nifti files
            digimouseAtlas_nii = load_untouch_nii(char(file(o,"DIGIMOUSE")));
            allenAtlas_nii = load_untouch_nii(char(file(o,"ALLEN")));

            % EF results
            efM_nii = load_untouch_nii(char(file(o,"EMAGNII")));
            efX_nii = load_untouch_nii(char(file(o,"EXNII")));
            efY_nii = load_untouch_nii(char(file(o,"EYNII")));
            efZ_nii = load_untouch_nii(char(file(o,"EZNII")));

            %% Find the voxels indices corresponding to digimouse atlas
            % We want to sample only from the digimouse brain.
            % The ef in voxel space (efm_nifti) has some extrapolated regions
            % that need to be disregarded.

            % We also split the digimouse atlas labels into left and right hemisphere.
            % (The original does not have this distinction).
            % This is necessary to analyze hemispheres individually given the interest
            % in the spatial effects of the stimulation montage
            % region id numbers
            % (TODO: make a cell for this info, and call it below)
            % medulla: left 4 --right 4+10
            % cerebellum: left 5 -- right 5+10
            % olfactory bulbs: left 6 --right 6+10
            % external cerebrum: left 7 -- right 7+10
            % striatum: left 8 -- right 8+10
            % rest of the brain: left 10 -- right 10+10
            %(right ids were made by adding the max id from the left, 10 here)

            idx_l = find(digimouseAtlas_nii.img>= 4 & digimouseAtlas_nii.img<= 10);
            idx_r = find(digimouseAtlas_nii.img>= 14 & digimouseAtlas_nii.img<= 20);
            switch (pv.hemisphere)
                case "left"
                    idx = idx_l;
                case "right"
                    idx = idx_r;
                case "both"
                    idx = [idx_l;idx_r];
            end
            brain_size = numel([idx_l;idx_r]);

            %% Mask Allen Atlas and EF matrices using digimouse Atlas voxels of interest
            % i.e., extract the voxels of interest for the other matrices:
            % for Allen Atlas (left, right or both)
            allenAtlas = allenAtlas_nii.img(idx);
            % for EF images (3D arrays)
            efM= efM_nii.img(idx);
            efX= efX_nii.img(idx);
            efY= efY_nii.img(idx);
            efZ= efZ_nii.img(idx);       
            inArea  = ismember(allenAtlas,allenIds);
            % extract the EF magnitude and components.
            efM= efM(inArea);
            efX= efX(inArea);
            efY= efY(inArea);
            efZ= efZ(inArea);
            
            % Determine summary
            S = EFMouse.summary([efX efY efZ efM]);
            brain_percentage = 100*sum(inArea)/brain_size;
            fprintf('Area: %s (%.1f%% of brain) , hemisphere %s \n',area,brain_percentage,pv.hemisphere)
            disp(S)
            
            % Determine focality
            % make a 3D mask for area of interest
            x = allenAtlas_nii.img;
            area_mask = zeros(size(x));
            area_mask(idx) = 1; % hemisphere mask: left, right or both            
            z = ismember(x,allenIds); % extract area of interest voxel indices (no hemispheric difference, ie. 10 is both left and right)
            area_mask = area_mask & z; % intersect hemispheric mask and z
            efX= efX_nii.img(area_mask);
            efY= efY_nii.img(area_mask);
            efZ= efZ_nii.img(area_mask);
            area_ef = [efX,efY,efZ];

            % make a 3D mask for the reference area
            ref_mask = zeros(size(x));
            ref_mask([idx_l;idx_r]) = 1; % get both hemispheres
            [~,referenceIds] = queryAllen(o,pv.foc_reference);
            z = ismember(x,referenceIds);
            ref_mask = ref_mask & z; % intersect hemispheric mask and z
            % find the set difference: use index of the logicals
            ref_ids = setdiff(find(ref_mask),find(area_mask));
            efX= efX_nii.img(ref_ids);
            efY= efY_nii.img(ref_ids);
            efZ= efZ_nii.img(ref_ids);
            ref_ef = [efX,efY,efZ];
            %%% for debugging, 3D mask of the reference area
            %ref_mask = zeros(size(x));
            %ref_mask(ref_ids) = 1;

            F = EFMouse.focality(area_ef,ref_ef,threshold = pv.foc_threshold,percentile_max = pv.foc_percentile_max);
            fprintf('Relative focality = %.4f%% of %d reference voxels (eMag > %.2f%% of the area max (%.2fth percentile))\n',F(1)*100,F(2),pv.foc_threshold,pv.foc_percentile_max);

            % Compute homogeneity
            H = EFMouse.homogeneity(area_ef);
            fprintf('[Homogeneity ranges from 0 to 1]\n');
            fprintf('Homogeneity = %.4f, for ef_norm_mean: (%.3f %.3f %.3f)\n',H(1),H(2),H(3),H(4));

            fprintf('----analyzeAtlas elapsed time: %.4f seconds\n',toc);

        end
        function [names,ids,acronyms,RGB] = queryAllen(o,area)
            % queryAllen retrieves information on an area of interest
            % in the Allen atlas.
            % area: A single name identifying an region in the Allen Atlas.
            % Output is a list of names, ids (to identify voxels in the
            % Allen Atlas Nifti file), and the corresponding acronyms and
            % colors.
            % aux_files/allenAtlas_labels.mat was created from a cell array
            % exported from the Allen Atlas website.
            arguments
                o (1,1) EFMouse
                area (1,1) string
            end

            % Read the hierarcy and labels from a file (constructed from
            % the Allen Atlas output files).
            load(file(o,"ALLENLABELS"),'hierarchy','labels');
            if ~ismember(area,labels.Name)
                fprintf("%s is not an area in the Allen Atlas.\n ",area);
                names= "";
            else
                % find the row and col of all areas belonging the area_of_interest
                [row,col] = find(strcmpi(hierarchy,area));

                % loop and cache every leave of the tree starting in area_of_interest
                % the data structure is organized in a tree form
                names = [];
                for r = 1:length(row)
                    elements = hierarchy(row(r),col(r)+1:end);
                    for e = 1:length(elements)
                        %if isempty(elements(e));continue;end
                        names = [names,string(elements{e})]; %#ok<AGROW>
                    end

                end
                names = unique(names)';
            end
            if nargout>1
                stay = ismember(labels.Name,names);
                ids = labels.ID(stay);
                if nargout >2
                    acronyms = labels.Acronym(stay);
                    if nargout >3
                        RGB = labels.RGB(stay);
                    end
                end
            end
        end


        function output_aligned = align2DigimouseAtlas_fsl(o,fNameUnaligned)
            % apply FSL linear transformation and trilinear interpolation
            % requires FSL installation
            % Uses a transformation matrix that aligns to the
            % aux_files/digimouseAtlas.nii.gz
            % calls FSL from terminal using matlab "system"
            arguments
                o (1,1) EFMouse
                fNameUnaligned (1,1) string
            end

            % Transformation matrix
            % (this was made using landmark based alignment with rotations,
            % and translations. The olfactory bulbs and the medulla served as
            % references.
            trans_matrix =file(o,"TRANS"');
            % Reference nifti file
            % FSL takes header information from this one, vox dim and othe important
            % parameters
            reference_nii = file(o,"DIGIMOUSE");

            % Define output file name
            % Delete the unAligned part of the name
            output_aligned = erase(fNameUnaligned,'_unaligned');


            % Find FSL directory
            % call from terminal
            if ispc
                % Using Windows Subsystem for Linux to run FSL
                % Assuming that the FSLDIR is ~/fsl (couldn't retrieve
                % FSLDIR from a command prompt).
                flirtCmd = 'wsl --shell-type login  ~/fsl/bin/flirt ';
                % Remap the filenames for WSL
                wslName = @(x) (strrep(strrep(x,extractBefore(x,3),"/mnt/" + lower(extractBefore(x,2))),"\","/"));
                trans_matrix = wslName(trans_matrix);
                reference_nii = wslName(reference_nii);
                output_aligned = wslName(output_aligned);
                fNameUnaligned = wslName(fNameUnaligned);
            else
                FSL_dir = getenv('FSLDIR');
                flirtCmd = [FSL_dir '/bin/flirt'];
            end

            % Build a string that will be the call for FSL
            % FSL needs full paths for all inputs and outputs!
            cmd= sprintf('%s -in "%s"  -applyxfm -init "%s" -out "%s" -paddingsize 0.0  -interp trilinear -ref "%s"', ...
                flirtCmd,fNameUnaligned, trans_matrix,output_aligned,reference_nii);
            % Run FSL call from terminal
            [status,cmdout] =system(cmd,'-echo');
            assert(status==0,"FLS call failed: %s",cmdout)
            fprintf('---%s created\n',output_aligned)
        end


    end

    

    %% Static helper functions
    methods (Static, Access=public)



        function nii = niftiCreator(input_3Dmatrix)
            % function to create nifti files specialized for our stimulation montage
            % modelling.
            % It will set the image such that:
            % voxel dim is 0.1 millimeters
            % and it is centered (0,0,0) in (medial,interaural,bregma)
            % following stereotaxic-Paxinos approach

            % make a nii file from ef_mag_F
            nii = make_nii(input_3Dmatrix);
            % Assign header parameters to the nifti file, such that it is
            % defined in the same coordinate space as aux_files/digimouseAtlas.nii.gz
            % these flags determine the "real world coordinates"
            nii.hdr.hist.sform_code = 2;
            nii.hdr.hist.qform_code = 0;
            % set the voxel dimensions to the digimouseAtlas 0.1 mm x voxel (iso)
            nii.hdr.dime.pixdim = [0,0.1000,0.1000,0.1000,0,0,0,0];
            %scaling = 1 avoids rescaling when displaying
            nii.hdr.dime.scl_slope = 1;
            % this flag defines the voxel units in  millimeters for the header
            nii.hdr.dime.xyzt_units = 2;
            % Use this to recenter the real world coordinates.
            % This is used to center the image according to our stereotaxic-Paxinos
            % approach
            % (0,0,0) real world --> (53,84,77) voxel space
            % (medial,interaural,bregma)
            nii.hdr.hist.srow_x = [0.1,0,0,-5.3];
            nii.hdr.hist.srow_y = [0,0.1,0,-8.4];
            nii.hdr.hist.srow_z = [0,0,0.1,-7.7];

        end

        function voxel_coordinates = mesh2VoxelCoordinates(X_target_mesh,Y_target_mesh,Z_target_mesh)
            % mesh2VoxelCoordinates transform mesh space coordinates to voxel space
            % coordinates. Inverse of voxel2meshCoordinates.m
            % Only use aux_files/EFMouse_digimouse_vox2mesh.nii.gz for the
            % voxel space coordinates!
            % We strongly recommend FSLeyes to open the nifti file and select the
            % target voxel coordinates.
            %
            % eg. voxel_coordinates = mesh2voxelCoordinates(-1.9290,31.6855,4.6887)

            x_min_voxel = 1;
            y_min_voxel = 1;
            z_min_voxel = 1;

            x_max_voxel = 380;
            y_max_voxel = 992;
            z_max_voxel = 208;


            x_min_mesh = -17.8017;
            y_min_mesh = -43.4963;
            z_min_mesh = -10.1639;

            x_max_mesh = 16.0450;
            y_max_mesh = 45.9763;
            z_max_mesh = 10.9903;


            % Perform linear transformation from mesh to voxel space
            X_target_voxel = round(((X_target_mesh - x_min_mesh) / (x_max_mesh - x_min_mesh)) * (x_max_voxel - x_min_voxel) + x_min_voxel);
            Y_target_voxel = round(((Y_target_mesh - y_min_mesh) / (y_max_mesh - y_min_mesh)) * (y_max_voxel - y_min_voxel) + y_min_voxel);
            Z_target_voxel = round(((Z_target_mesh - z_min_mesh) / (z_max_mesh - z_min_mesh)) * (z_max_voxel - z_min_voxel) + z_min_voxel);

            voxel_coordinates = [X_target_voxel,Y_target_voxel,Z_target_voxel];
        end

        function mesh_coordinates = voxel2MeshCoordinates(X_target_voxel,Y_target_voxel,Z_target_voxel)
            % voxel2MeshCoordinates transform voxel space coordinates to mesh space coordinates
            % to use for electrodes and craniotomy positioning or to define regions of interest
            % for analysis.
            %
            % Only use aux_files/EFMouse_digimouse_vox2mesh.nii.gz for the
            % voxel space coordinates!
            % We strongly recommend FSLeyes to open the nifti file and select the
            % target voxel coordinates.
            %
            % eg. mesh_coordinates = voxel2meshCoordinates(163,753,155)
            %

            %%
            x_min_voxel = 1;
            y_min_voxel = 1;
            z_min_voxel = 1;

            x_max_voxel = 380;
            y_max_voxel = 992;
            z_max_voxel = 208;


            x_min_mesh = -17.8017;
            y_min_mesh = -43.4963;
            z_min_mesh = -10.1639;

            x_max_mesh = 16.0450;
            y_max_mesh = 45.9763;
            z_max_mesh = 10.9903;

            % Perform linear transformation from voxel space to mesh space
            X_target_mesh = ((X_target_voxel - x_min_voxel) / (x_max_voxel - x_min_voxel)) * (x_max_mesh - x_min_mesh) + x_min_mesh;
            Y_target_mesh = ((Y_target_voxel - y_min_voxel) / (y_max_voxel - y_min_voxel)) * (y_max_mesh - y_min_mesh) + y_min_mesh;
            Z_target_mesh = ((Z_target_voxel - z_min_voxel) / (z_max_voxel - z_min_voxel)) * (z_max_mesh - z_min_mesh) + z_min_mesh;

            mesh_coordinates = [X_target_mesh,Y_target_mesh,Z_target_mesh];

        end

        function S = summary(ef)
            arguments
                ef (:,:) double  % Electric field estimates  [N 3] (x,y,z)  or [N 4] (x,y,z, magnitude)
            end
            if size(ef,2)==3
                % Compute magnitude from components
                ef = [ef sqrt(sum(ef.^2,2))];
            end
            components = num2cell(ef,1);
            varNames = {'eX','eY','eZ','eMag'};            
            T = table(components{:},'VariableNames',varNames);
            summaryFun = {@mean,@median,@std,@(x) min(nonzeros(x)),@max};
            summaryName = {'mean','median','std','min','max'};
            S =  table('RowNames',varNames);
            for i = 1:numel(summaryFun)
                tmp = varfun(summaryFun{i},T);
                S = addvars(S,tmp.Variables','NewVariableNames',summaryName{i});
            end
        end

        function F = focality(area_ef,reference_ef,pv)
            % A = number of elements in reference > x% of max of roi
            % B = total number of elements in complement
            % Focality = A/B
            % Lower number implies greater focality
            % Adapted from Fernandes et al. https://dx.doi.org/10.1088/1361-6560/ad222d
            % the defaults are from Fernandes et al.
            arguments
                area_ef (:,:) double
                reference_ef (:,:) double
                pv.threshold = 75 % 75% of the max
                pv.percentile_max = 99.9 % max defined as the 99.9 percentile
            end
            % compute x% max of roi as reference (as in Fernandes et al.)
            area_efMag = sqrt(sum(area_ef.^2,2));
            max_limit = prctile(area_efMag,pv.percentile_max); 
            cut_off = (pv.threshold/100) * max_limit;
            % now use the cutoff to threshold the complement_ef.
            % this is all the cortical area outside the roi
            reference_efMag = sqrt(sum(reference_ef.^2,2));
            num_nodes_p = sum(reference_efMag > cut_off);
            % compute focality
            reference_num_nodes = numel(reference_efMag);
            F = num_nodes_p/reference_num_nodes;
            % include focality and number of total nodes for reference
            F = [F,reference_num_nodes];
        end

        function H = homogeneity(ef)
            % Measure of homogeneity of field direction
            % goes from 0 to 1, from low to high homogeneity
            ef_mag = sqrt(sum(ef.^2,2));
            ef_norm = ef./ef_mag;
            ef_norm_mean = mean(ef_norm,1,"omitnan");
            H = sqrt(sum(ef_norm_mean.^2,2));
            H = [H,ef_norm_mean];
        end

    end
end