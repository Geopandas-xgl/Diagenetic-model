function dy = DiagModel_GLX_dJ6dt(t,J,m,v,F_init)
%% ==================== Initialization =====================
box = length(m);         % number of boxes
dy = zeros(v*box,1);     % Create a vector which will store all variables





fluid(1).C = F_init(1);  % initialze fluid chemistry vectors
fluid(1).O = F_init(2);
fluid(1).dC = F_init(3);
fluid(1).dO = F_init(4);

%% ============= set up a system of differential equations  ===============
for  n = 1:box

    i = v*n - v; 

    Prim.C = J(i+1);
    Prim.O = J(i+2);
    Secm.C = J(i+3);
    Secm.O = J(i+4);
    Prim.dC = J(i+5);
    Prim.dO = J(i+6);
    Secm.dC = J(i+7)/Secm.C; 
    Secm.dO = J(i+8)/Secm.O;
    fluid(n+1).C = J(i+11);
    fluid(n+1).O = J(i+12);
    fluid(n+1).dC = J(i+13)/fluid(n+1).C;
    fluid(n+1).dO = J(i+14)/fluid(n+1).O;
    
    %%%%%%% Set up function that solves all input and output fluxes to each box
    [in,out] = DiagModel_GLX_fluxes(Prim, fluid, n, box);
    
    %%%%%%% dM/dt = input fluxes - output fluxes
    dy(i+1) = in.Prim.C - out.Prim.C;
    dy(i+2) = in.Prim.O - out.Prim.O;
    dy(i+3) = in.Secm.C - out.Secm.C;
    dy(i+4) = in.Secm.O - out.Secm.O;
    dy(i+5) = 0;
    dy(i+6) = 0;
    dy(i+7) = in.Secm.dC - out.Secm.dC;
    dy(i+8) = in.Secm.dO - out.Secm.dO;
    dy(i+9) = in.solid.dC - out.solid.dC;
    dy(i+10) = in.solid.dO - out.solid.dO;
    dy(i+11) = in.fluid.C - out.fluid.C;
    dy(i+12) = in.fluid.O - out.fluid.O;
    dy(i+13) = in.fluid.dC - out.fluid.dC;
    dy(i+14) = in.fluid.dO - out.fluid.dO;
end
end
