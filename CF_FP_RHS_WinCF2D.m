% function to define the RHS of the ODEs 
% U' = F(U,SO,SN, F, P), here U = (theta_f; theta_p) 
% array of size (Nx+1) by 2


function Y = CF_FP_RHS_WinCF2D(U,SO, SN, SA,SG,F,P,I,TB,TT)

global mu_f mu_pa mu_pn K_f K_p K_O K_N  K_A K_G  beta0 beta1 beta2 C_0 d_p d_T Ind_T Ind_B

g0 = @(x) mu_f*(1-(2/pi)*atan(beta0*x));
g1 = @(x) K_f*(1 - 0.9*(2/pi)*atan(beta1*x));
g2 = @(x) K_p*(0.9*(0.5 - (1/pi)*atan(beta2*x)) + 0.1);

Y = zeros(size(U));

Y(:,1) = g0(SO).*(SG./(K_G + SG)).*U(:,1).*(1-U(:,1)./g1(I));  

Y(:,2) = (mu_pa*(SO./(K_O + SO)) + exp(-d_p*SO)*mu_pn.*(SN./(K_N + SN))).*(SA./(K_A + SA)).*U(:,2).*(1 - U(:,2)./g2(F - P + C_0 - Ind_B*TB)) - Ind_T*d_T*TT.*U(:,2);

end





