classdef Stage < double
    % Enumeration class used to order operations in EFMouse.run
    enumeration
        UNDEF (-1) % Undefined -  newly created objects start here
        INIT (0)   % Initialized - base model setup
        MESH (1)   % Meshed - mesh modified to include electrodes/craniotomies
        EXPORT (2) % Export - Save .pro and .msh files for getdp
        GETDP (3)  % Run GetDP
        ATLAS (4)  % Save results as .nii.gz files for use in atlas based analysis.
    end

end