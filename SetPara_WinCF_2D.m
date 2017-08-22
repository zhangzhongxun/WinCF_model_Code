%parameter values for 2D WinCF model

global K_f K_p K_O K_N K_G K_A D_F D_P D_O D_N D_I D_A D_G D_B D_T Y_pPo Y_pPn  Y_pO Y_pN Y_pIo Y_pIn Y_pA  Y_fF Y_fG d_O d_p d_T Y_pT % dimensionless parameters, global variables

global mu_f mu_pa mu_pn  %growth rates of bacteria

global NX NY xL xR yB yT dt NT FS_interval t_Final % spatial and temporal domain and mesh size 

global beta0 beta1 beta2 C_0 flag_BC Ind_B Ind_T % control parameters 
 
global SO_BC0  theta_f_0 theta_p_0 F_IC0 P_IC0 I_IC0  SN_IC0 SA_IC0 SG_IC0 TB_IC0 TT_IC0


t0 = 3.6e3; %reference time scale, 1 hour
h0 = 1e-2; % reference length scale, 100 micro-meter

mu_f0 = .2/3600; % .2/3600 %dimensional growth rate of fermenter
K_f0 = 2; %dimensional maximum carrying capacity of fermenter
theta_f0 = 1; % dimensional fermenter cell density

mu_pa0 = .9/3600; %dimensional growth rate of Pseudomonas by aerobic respiration
mu_pn0 = .3/3600; %dimensional growth rate of Pseudomonas by anaerobic respiration
d_p0 = 2; % dimensional inhibition rate of Pseudomonas anaerobic respiration by oxygen
K_p0 = 5; %dimensional maximum carrying capacity of Pseudomonas
theta_p0 = 1; % dimensional Pseudomonas cell density
d_T0 = 1.8e-4; % dimensional killing rate of Pseudomona

S_O0 = 1; % reference oxygen concentration
S_N0 = 1; % reference nitrate concentration
S_A0 = 1; % reference amino acid concentration
S_G0 = 1; % reference sugar concentration
P0 = 1; % reference ammonium concentration
F0 = 1; % reference acid concentration
I0 = 1; % reference fermenter inhibitor concentration
TB0 = 1; % reference bicarbonate concentration
TT0 = 1; % reference tobramycin concentration

K_O0 = .5; % dimensional half-saturation constant for Monod growth kinetic of Pseudomonas on oxygen
K_N0 = .5;  % dimensional half-saturation constant for Monod growth kinetic of Pseudomonas by anaerobic respirationpp
K_A0 = .8; % dimensional half-saturation constant for Monod growth kinetic of Pseudomonas on amino acid
K_G0 = .9; % dimensional half-saturation constant for Monod growth kinetic of Fermenters on sugar


D_F0 = 3e-10; %dimensional diffusion coefficient for acid

D_P0 = 1e-9; %2e-10; %dimensional diffusion coefficient for ammonium

D_O0 = 1.97e-9; %dimensional diffusion coefficient for oxygen

D_N0 = 1e-10; %dimensional diffusion coefficient for nitrate

D_I0 = 3e-10; %1e-10; %dimensional diffusion coefficient for fermenter inhibitor

D_A0 = 8e-10; %dimensional diffusion coefficient for amino acid

D_G0 = 5.7e-10; %dimensional diffusion coefficient for sugar

D_B0 = 1.4e-10; %dimensional diffusion coefficient for bicarbonate

D_T0 = 1.5e-10; % 3.9e-10; %dimensional diffusion coefficient for tobramycin


%r_F0 = 2; % 1; %dimensional production rate of acid by fermenter

Y_pPo0 = 0.03; % dimensional yield coeff of pseudomonas on ammonium on aerobic respiration

Y_pPn0 = 0.05; % dimensional yield coeff of pseudomonas on ammonium on anaerobic respiration

Y_pO0 = 0.015; % dimensional yield coeff of pseudomonas on oxygen

Y_pN0 = 0.38; % dimensional yield coeff of pseudomonas on nitrate

Y_pIo0 = 0.3;% 0.5; % dimensional yield coeff of pseudomonas on fermenter inhibitor on aerobic respiration

Y_pIn0 = 3; %0.35; % dimensional yield coeff of pseudomonas on fermenter inhibitor on anaerobic respiration

Y_pA0 = 0.38; % dimensional yield coeff of pseudomonas on amino acid

Y_fF0 = .02; % dimensional yield coeff of fermenter on acid

Y_fG0 = 0.5; % 0.22; % dimensional yield coeff of fermenter on sugar

Y_pT0 = 0.5; % dimensional yield coeff of pseudomonas on tobramycin

d_O0 = 1.5/3600; % dimensional oxygen consumption rate for endogenous respiration of Pseudomonas

beta0_0 = 10; % dimensional scaling coefficient for inhibition of fermenter by oxygen
beta1_0 = 10; % 10; % dimensional scaling coefficient for inhibition of fermenter by I, note beta1 = 0 means I has no effect
beta2_0 = 1;  % 1 % dimensional scaling coefficient for inhibition of Pseudomonas by lower pH

%%%%%%%%%%%%%%% Dimensionless parameter values  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_f = mu_f0*t0;     %dimensionless growth rate of fermenter
K_f = K_f0/theta_f0; %dimensionless  maximum carrying capacity of fermenter 
K_G = K_G0/S_G0;     % dimensionless half saturation constant of fermenter on sugar
beta0 = beta0_0*S_O0;
beta1 = beta1_0*I0;

mu_pa = mu_pa0*t0;  %dimensionless aerobic growth rate of Pseudomonas 
mu_pn = mu_pn0*t0;  %dimensionless anaerobic growth rate of Pseudomonas 
d_p = d_p0*S_O0; % dimensionless inhibition rate of Pseudomonas anaerobic respiration by oxygen
K_O = K_O0/S_O0;  % dimensionless half-saturation constant for Monod growth kinetic of Pseudomonas on oxygen
K_N = K_N0/S_N0; % dimensionless half-saturation constant for Monod growth kinetic of Pseudomonas by anaerobic respiration
K_A = K_A0/S_A0;  % dimensionless half-saturation constant for Monod growth kinetic of Pseudomonas on amino acid
K_p = K_p0/theta_p0; %dimensionless  maximum carrying capacity of Pseudomonas
beta2 = beta2_0*F0;
d_T = d_T0*t0*TT0; % dimensional killing rate of Pseudomona

D_F = D_F0*t0/h0/h0; %dimensionless diffusion coefficient for acid
Y_fF = Y_fF0*F0/theta_f0; % dimensionless yield coeff of fermenter on acid

D_P = D_P0*t0/h0/h0;%dimensionless diffusion coefficient for ammonium
Y_pPo = Y_pPo0*P0/theta_p0; % dimensionless yield coeff of pseudomonas on ammonium on aerobic respiration
Y_pPn = Y_pPn0*P0/theta_p0; % dimensionless yield coeff of pseudomonas on ammonium on anaerobic respiration


D_O = D_O0*t0/h0/h0; %dimensionlessl diffusion coefficient for oxygen
Y_pO = Y_pO0*S_O0/theta_p0; % dimensionless yield coeff of pseudomonas on oxygen


D_N = D_N0*t0/h0/h0; %dimensionless diffusion coefficient for nitrate
Y_pN = Y_pN0*S_N0/theta_p0; % dimensionless yield coeff of pseudomonas on nitrate

D_I = D_I0*t0/h0/h0; %dimensionless diffusion coefficient for fermenter inhibitor
Y_pIo = Y_pIo0*I0/theta_p0;  % dimensionless yield coeff of pseudomonas on fermenter inhibitor on aerobic respiration
Y_pIn = Y_pIn0*I0/theta_p0;  % dimensionless yield coeff of pseudomonas on fermenter inhibitor on anaerobic respiration

D_A = D_A0*t0/h0/h0; %dimensionless diffusion coefficient for amino acid
Y_pA = Y_pA0*S_A0/theta_p0; % dimensionless yield coeff of pseudomonas on amino acid

D_G = D_G0*t0/h0/h0; %dimensionless diffusion coefficient for sugar
Y_fG = Y_fG0*S_G0/theta_f0; % dimensionless yield coeff of fermenter on sugar

d_O = d_O0*t0; %  dimensionless oxygen consumption rate for endogenous respiration of Pseudomonas

C_0 = 5; % 0: high pH, 5: medium pH, 10: low pH; % background  value for pH, larger C_0 means lower pH

D_B = D_B0*t0/h0/h0; %dimensionless diffusion coefficient for bicarbonate
D_T = D_T0*t0/h0/h0; %dimensionless diffusion coefficient for tobramycin
Y_pT = Y_pT0*TT0/theta_p0; %dimensionless yield coeff of pseudomonas on tobramycin

%%%%%%%%%%%%% Dirichlet BC values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_BC = 2; % 1: in vivo; 2: WinCF tube

SO_BC0 = 1;

theta_f_0 = 0.1; 
theta_p_0 = 0.2;

F_IC0 = 0;
P_IC0 = 0; 
I_IC0 = 0;  
SA_IC0 = 1;
SG_IC0 = 1;
SN_IC0 = 1;

TB_IC0 = 2; % total amount of bicarbonate added to the tube
TT_IC0 = 2e-2; % total amount of tobramycin added to the tube

Ind_B = 1;  % Indicator for bicarbonate treatment
Ind_T = 1;  % Indicator for bicarbonate treatment

%%%%%%%%%%%%%  domain information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL = 0;
xR = .1;
yB = 0;
yT = .4;

NX = 5; % 20;
NY = 20; % 80;

t_Final = 50;  %  total simulation simulation time

dt = 0.05; % time step size

NT = floor(t_Final/dt);

Nstride = 20; % 10; %  Number of saved steps

FS_interval = floor(NT/Nstride); % File saving interval







