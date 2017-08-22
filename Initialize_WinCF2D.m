%function to initialize the unknowns
function [theta_f_n, theta_p_n, F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, TB_n, TT_n] = Initialize_WinCF2D(NNM)

global theta_f_0 theta_p_0

global F_IC0 P_IC0 I_IC0 SO_BC0 SN_IC0 SA_IC0 SG_IC0 TB_IC0 TT_IC0

global NX NY xL xR yB yT

dy = (yT - yB)/NY; % mesh size in y-direction

theta_f_n = theta_f_0*ones(NNM,1);

theta_p_n = theta_p_0*ones(NNM,1);

F_n = F_IC0*ones(NNM,1);

P_n = P_IC0*ones(NNM,1);

I_n = I_IC0*ones(NNM,1);

SO_n = SO_BC0*ones(NNM,1);

SN_n = SN_IC0*ones(NNM,1);

SA_n = SA_IC0*ones(NNM,1);

SG_n = SG_IC0*ones(NNM,1);

TB_n = (TB_IC0/(dy*(xR - xL)))*ones(NNM,1);
TB_n(1:(NX+1)*(NY-1)) = 0;

TT_n = (TT_IC0/(dy*(xR - xL)))*ones(NNM,1);
TT_n(1:(NX+1)*(NY-1)) = 0;

end

 
 