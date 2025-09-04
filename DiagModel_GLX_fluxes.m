function [in, out] = DiagModel_GLX_fluxes(Prim, fluid, n, box)

%%%%%%% Set up global structures
global alpha 
global R 
global Box 
global V 

%%%%%%% Set up linear variation of reaction rate (R) and floe rate (V)
R.Prim = -(R.start-R.end)/(box)*n + (R.start+(R.start-R.end)/(box)); 
V.equi = (-(V.start-V.end)/(box)*n + (V.start+(V.start-V.end)/(box)))/Box.length;

%%%%%% Flux into diagenetic mineral
dia.Prim.C = R.Prim*Prim.C;
dia.Secm.C = dia.Prim.C;

dia.Prim.O = R.Prim*Prim.O;
dia.Secm.O = dia.Prim.O;

%%%%%% Sediment input and output fluxes (F_in - F_out)
in.Prim.C = 0;
out.Prim.C =dia.Prim.C;
in.Prim.O = 0;
out.Prim.O =dia.Prim.O;

in.Secm.C = dia.Secm.C;
out.Secm.C =0;
in.Secm.O = dia.Secm.O;
out.Secm.O =0;

dia.c.C = fluid(n+1).dC + 1e3*log(alpha.C);
dia.c.O = fluid(n+1).dO + 1e3*log(alpha.O);

in.Secm.dC = in.Secm.C*dia.c.C;
out.Secm.dC = 0;
in.Secm.dO = in.Secm.O*dia.c.O;
out.Secm.dO = 0;

in.solid.dC = in.Secm.C*dia.c.C;
out.solid.dC = out.Prim.C*Prim.dC;
in.solid.dO = in.Secm.O*dia.c.O;
out.solid.dO = out.Prim.O*Prim.dO;

in.fluid.C = V.equi*fluid(n).C + out.Prim.C;
out.fluid.C = V.equi*fluid(n+1).C + in.Secm.C;
in.fluid.O = V.equi*fluid(n).O + out.Prim.O;
out.fluid.O = V.equi*fluid(n+1).O + in.Secm.O;


in.fluid.dC = V.equi*fluid(n).C*fluid(n).dC + out.solid.dC;
out.fluid.dC = V.equi*fluid(n+1).C*fluid(n+1).dC + in.Secm.dC;

in.fluid.dO = V.equi*fluid(n).O*fluid(n).dO + out.solid.dO;
out.fluid.dO = V.equi*fluid(n+1).O*fluid(n+1).dO + in.Secm.dO;
end
