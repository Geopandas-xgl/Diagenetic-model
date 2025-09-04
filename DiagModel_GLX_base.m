%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   DiagModel_Initial function    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Solid,WR] = DiagModel_GLX_base(index)
%% ==================== Define parameters =====================

%%%%%%% Set up global structures
global Fluid 
global alpha 
global R 
global Box  
global box 
global V

%%%%%%% Set up baseline conditions
PrimC = [3,-4.5];               % Primary calcite dC (permil) (Normalized to V-PDB)  
PrimO = [0.5,-1];               % Primary calcite dO (permil) (Normalized to V-PDB) 
dC = [-14, -14];                % d13C_fluid values (permil)
dO = [-46,-46];                 % d18O_fluid_values (permil)
SecMineral = 0.001;             % Secondary calcite fraction in original carbonate
PriMineral = 1-SecMineral;      % Primary calcite fraction in original carbonate 
Solid.Prim.dC = PrimC(index);   % Prim calcite dC (permil) 
Solid.Prim.dO = PrimO(index);   % Prim calcite dO (permil)
C_fluid = 0.0012;               % [C] in fluid (mol/kg) 
O_fluid = 1000/18;              % [O] in fluid (mol/kg)
Molar.C = 12.011;               % C molar mass (g/mol)
Molar.O = 15.999;               % O molar mass (g/mol)
Molar.Prim = 100.086;           % Primary calcite molar mass (g/mol)
Molar.Secm = 100.086;           % Second calcite molar mass (g/mol)
Conc.Prim.C = Molar.C/Molar.Prim*1e6;   % [C] in primary calcite (ug/g)
Conc.Prim.O = Molar.O*3/Molar.Prim*1e6; % [O] in primary calcite (ug/g)
Conc.Secm.C = Molar.C/Molar.Secm*1e6;   % [C] in secondary calcite (ug/g)
Conc.Secm.O = Molar.O*3/Molar.Secm*1e6; % [O] in secondary calcite (ug/g)
Fluid.dC = dC(index);                   % d13C_fluid values (permil)
Fluid.dO = dO(index);                   % d18O_fluid_values (permil)
alpha.C = 0.9988;                       % d13C_alpha; Isotopic fractors of diageneetic mineral
alpha.O = 1.031;                        % d18C_alpha; Isotopic fractors of diageneetic mineral  
Box.length = 1;                         % length of a box (m)                 
Box.vol = (Box.length*100)^3;           % Volume of a box (cm^3)
Box.Phi = 0.3;                          % Pore volume fraction (or Porosity)  
Box.RhoS = 1.8;                         % Density of sediment (g/cm3) 
Box.RhoF = 1.0;                         % Density of fluid (g/cm3)
Box.sed = (1-Box.Phi)*Box.vol*Box.RhoS/1e3;       % Sediment mass (kg)
Box.fluid = Box.Phi*Box.vol*Box.RhoF/1e3;         % Fluid mass (kg) 
V.start = 1;                                      % Flow Rate for the first box (m/yr-1) 
V.end = 0.1;                                      % Flow Rate for the end box (m/yr-1)
R.start = 1e-5;                                   % Reaction rate for the first box (g/g yr-1)
R.end = 0.2e-5;                                   % Reaction rate for the end box (g/g yr-1)
time = 0.5e6;                                     % Reaction time (yr)
box = 50;                                         % Box number


Solid.Prim.C = PriMineral*Box.sed*Conc.Prim.C*1e-3/Molar.C; % moels C in parimary calcite (mol)
Solid.Prim.O = PriMineral*Box.sed*Conc.Prim.O*1e-3/Molar.O; % moels O in parimary calcite (mol)
Solid.Secm.C = SecMineral*Box.sed*Conc.Secm.C*1e-3/Molar.C; % moels C in secondary calcite (mol)
Solid.Secm.O = SecMineral*Box.sed*Conc.Secm.O*1e-3/Molar.O; % moels O in secondary calcite (mol)
Solid.Secm.dC = (Fluid.dC + 1e3*log(alpha.C))*Solid.Secm.C; % d13C in secondary calcite (permil)
Solid.Secm.dO = (Fluid.dO + 1e3*log(alpha.O))*Solid.Secm.O; % d18O in secondary calcite (permil)
Solid.dC = (Solid.Prim.C*Solid.Prim.dC + Solid.Secm.dC);    % d13C in bulk carbonate (permil)
Solid.dO = (Solid.Prim.O*Solid.Prim.dO + Solid.Secm.dO);    % d18O in bulk carbonate (permil)
Fluid.C = C_fluid*Box.fluid;                                % moles C in fluid (mol)
Fluid.O = O_fluid*Box.fluid;                                % moles O in fluid (mol)
Fluid.dC_init = Fluid.C*Fluid.dC;                           % Initial (moles C)*d13C of fluid to solve ODE 
Fluid.dO_init = Fluid.O*Fluid.dO;                           % Initial (moles O)*d18O of fluid to solve ODE 

%% ==================== Initial conditionas vetor for ODE =====================
par = [Solid.Prim.C; Solid.Prim.O; Solid.Secm.C; Solid.Secm.O; 
       Solid.Prim.dC; Solid.Prim.dO; Solid.Secm.dC; Solid.Secm.dO; 
       Solid.dC; Solid.dO; 
       Fluid.C; Fluid.O; Fluid.dC_init; Fluid.dO_init]; 

v = length(par);                               % Number of variables 
tspan = [0, time];                             % Time span (yr) 
m = 1 : box;                                   % Box vector  
F.V = [Fluid.C, Fluid.O,  Fluid.dC, Fluid.dO]; % Intial fluid vector
  
for i = 1:length(m)
    Jt0(:,i) = par(:,1);                       % Initialize all boxes 
end
%%  ==================== Evolve the model ===================== 

tic, [t,J] = ode15s(@(t,J)DiagModel_GLX_dJ6dt(t,J,m,v,F.V), tspan, Jt0); toc

%%  ================== Model ouput =============================
for i = 1:length(m) 
    j = v*i-v; 
    a = j + 1;
    c = j + v;
    Jprime(:,:,i) = J(1:length(t), a:c)'; % transpose J (time vs tracked variables) for box  i and
end
%%  ============== Extract variables from Jprime (3d matrix)=============
for i = 1:length(t) 
    for j = 1:length(m)
Solid.dC(i,j) = Jprime(9,i,j)/(Jprime(1,i,j)+Jprime(3,i,j));  % The d13C in bulk diagenetic carbonate (permil)
Solid.dO(i,j) = Jprime(10,i,j)/(Jprime(2,i,j)+Jprime(4,i,j)); % The d18O in bulk diagenetic carbonate (permil)
WR(i,j) = t(i);                                               % Diagenetic duration (yr)
    end
end
end