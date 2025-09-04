%%                          Carbonate Dissolution-Precipitation Diagenesis Model
%
%           JGM_V5.0: Dissolution of carbonate and precipitation of diagenetic mineral
%
%           2015-08-04      Initial Build                   Anne-Sofie C. Ahm (Princeton Univ)
%           2015-11-30      Revised to include Sr           Anne-Sofie C. Ahm (Princeton Univ)
%           2018-03-09      Published                       Ahm et al. (2018) GCA 236 (2018) 140-159
%           2018-11-07      Revised to include Li           Jack Geary Murphy  (Princeton Univ)
%           2020-07-15      Revised to include HMC        Jack Geary Murphy  (Princeton Univ)
%           2021-02-06      Revised to include prim Dolo  Jack Geary Murphy  (Princeton Univ)
%           2025-9-02       Revised to phreatic meteoric diagenesis  Guolin Xiong  (Nanjing university)
% =============================================================
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Run diagenetic model          %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clc

ModeData = struct("Solid",cell(1,2),"WR",cell(1,2));

parfor N = 1:2  % "N = 1" denotes inital dC and dO of primary minerals in Pre-CIE;
                % "N = 2" denotes inital dC and dO of primary minerals in Peak-CIE;
    [Solid,WR] = DiagModel_GLX_base(N); % Acquire diagenetic reults
    ModelData(N).Solid = Solid; 
    ModelData(N).WR = WR; 
end

% Save results
dataFolder = fullfile('Data');
if ~isfolder(dataFolder)
    mkdir(dataFolder);
end

save("Data\ModelData.mat","ModelData");
%% Plot 
DiagModel_GLX_plot

