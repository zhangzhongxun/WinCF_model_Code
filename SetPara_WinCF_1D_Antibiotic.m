%parameter values for 2D WinCF model

global K_f K_p K_O K_N K_G K_A D_F D_P D_O D_N D_I D_A D_G Y_pPo Y_pPn  Y_pO Y_pN Y_pIo Y_pIn Y_pA  Y_fF Y_fG d_O d_p 
global D_Tp D_Tf D_Tw d_Tp d_Tf d_pw d_fw Y_pT Y_fT Y_pw Y_fw % dimensionless parameters, global variables

global mu_f mu_pa mu_pn  %growth rates of bacteria

global L Nx dx dt NT FS_interval t_Final % spatial and temporal domain and mesh size 

global beta0 beta1 beta2 C_0 Ind_Tp Ind_Tf Ind_Tw  % control parameters 
 
global SO_BC1  theta_f_0 theta_p_0 F_IC0 P_IC0 I_IC0  SN_IC0 SA_IC0 SG_IC0 Tf_IC0 Tp_IC0 Tw_IC0

global BC_flag_F BC_flag_P BC_flag_I BC_flag_SO BC_flag_SN BC_flag_SA BC_flag_SG BC_flag_Tf BC_flag_Tp BC_flag_Tw

global Treat_ID  N_Treat Treat_Names Ratio_pf_ini N_ratio Sample_ID  % global variables for treatment combination and initial pseudomonas and fermenter ratio

global Save_Inter_flag

global lambda 

t0 = 3.6e3; %reference time scale, 1 hour
h0 = 1e-2; % reference length scale, 1 center meter

mu_f0 = .2/3600; % .2/3600 %dimensional growth rate of fermenter
K_f0 = 5; %dimensional maximum carrying capacity of fermenter ****** parameter changed ******
theta_f0 = 1; % dimensional fermenter cell density

mu_pa0 = .9/3600; %dimensional growth rate of Pseudomonas by aerobic respiration
mu_pn0 = .3/3600; %dimensional growth rate of Pseudomonas by anaerobic respiration
d_p0 = 2; % dimensional inhibition rate of Pseudomonas anaerobic respiration by oxygen
K_p0 = 5; %dimensional maximum carrying capacity of Pseudomonas
theta_p0 = 1; % dimensional Pseudomonas cell density
d_Tp0 = 8e-2; % 5e-3; % 1.8e-4; % dimensional killing rate of antibiotic only kills Pseudomona ****** parameter changed ******
d_Tf0 = 5e-3; % 1.8e-4; % dimensional killing rate of antibiotic only kills anarobes
d_pw0 = 5e-3; % 1.8e-4; % dimensional killing rate on pseudomonas of wide spectrum antibiotic
d_fw0 = 5e-3; % 1.8e-4; % dimensional killing rate on anarobes of wide spectrum antibiotic

S_O0 = 1; % reference oxygen concentration
S_N0 = 1; % reference nitrate concentration
S_A0 = 1; % reference amino acid concentration
S_G0 = 1; % reference sugar concentration
P0 = 1; % reference ammonium concentration
F0 = 1; % reference acid concentration
I0 = 1; % reference fermenter inhibitor concentration
Tp0 = 1; % reference antibiotic concentration only kills pseudomonas
Tf0 = 1; % reference antibiotic concentration only kills anarobes
Tw0 = 1; % reference antibiotic concentration wide spectrum

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

D_Tp0 = 1.4e-10; %dimensional diffusion coefficient for antibiotic only kills pseudomonas

D_Tf0 = 1.5e-10; %dimensional diffusion coefficient for antibiotic only kills anarobes

D_Tw0 = 1.5e-10; %dimensional diffusion coefficient for antibiotic kills both


%r_F0 = 2; % 1; %dimensional production rate of acid by fermenter

Y_pPo0 = 0.03; % dimensional yield coeff of pseudomonas on ammonium on aerobic respiration

Y_pPn0 = 0.05; % dimensional yield coeff of pseudomonas on ammonium on anaerobic respiration

Y_pO0 = 0.005; % 0.015; % dimensional yield coeff of pseudomonas on oxygen ****** parameter changed ******

Y_pN0 = 0.38; % dimensional yield coeff of pseudomonas on nitrate

Y_pIo0 = 0.3;% 0.5; % dimensional yield coeff of pseudomonas on fermenter inhibitor on aerobic respiration

Y_pIn0 = 3; %0.35; % dimensional yield coeff of pseudomonas on fermenter inhibitor on anaerobic respiration

Y_pA0 = 0.38; % dimensional yield coeff of pseudomonas on amino acid

Y_fF0 = .02; % dimensional yield coeff of fermenter on acid 

Y_fG0 = 2; % 0.5; % 0.22; % dimensional yield coeff of fermenter on sugar ****** parameter changed ******

Y_pT0 = 50; % 0.5; % dimensional yield coeff of pseudomonas on antibiotic only kills pseudomonas

Y_fT0 = 50; % 0.5;  % dimensional yield coeff of anarobes on antibiotic only kills anarobes

Y_pw0 = 50; % 0.5;  % dimensional yield coeff of pseudomonas on wide spectrum antibiotic

Y_fw0 = 50; % 0.5;  % dimensional yield coeff of anarobes on wide spectrum antibiotic

d_O0 = 1.5/3600; % dimensional oxygen consumption rate for endogenous respiration of Pseudomonas

beta0_0 = 10; % 10; % dimensional scaling coefficient for inhibition of fermenter by oxygen ****** parameter changed ******
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

D_Tp = D_Tp0*t0/h0/h0; %dimensionless diffusion coefficient for antibiotic only kills Pseudomona
D_Tf = D_Tf0*t0/h0/h0; %dimensionless diffusion coefficient for antibiotic only kills anarobes
D_Tw = D_Tw0*t0/h0/h0; %dimensionless diffusion coefficient for wide spectrum antibiotic
Y_pT = Y_pT0*Tp0/theta_p0; %dimensionless yield coeff of pseudomonas on antibiotic only kills Pseudomona
Y_fT = Y_fT0*Tf0/theta_f0; %dimensionless yield coeff of anarobes on antibiotic only kills anarobes
Y_pw = Y_pw0*Tw0/theta_p0; %dimensionless yield coeff of pseudomonas on wide spectrum antibiotic
Y_fw = Y_fw0*Tw0/theta_f0; %dimensionless yield coeff of anarobes on wide spectrum antibiotic


d_Tp = d_Tp0*t0*Tp0; % dimensionless killing rate of antibiotic only kills Pseudomona
d_Tf = d_Tf0*t0*Tf0; % dimensionless killing rate of antibiotic only kills anarobes
d_pw = d_pw0*t0*Tw0; % dimensionless killing rate on pseudomonas of wide spectrum antibiotic
d_fw = d_fw0*t0*Tw0; % dimensionless killing rate on anarobes of wide spectrum antibiotic


%%%%%%%%%%%%% Dirichlet BC values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BC_flag_F = 1;
BC_flag_P = 1;
BC_flag_I = 1; 
BC_flag_SO = 2;
BC_flag_SN = 1;
BC_flag_SA = 1; 
BC_flag_SG = 1;
BC_flag_Tf = 1;
BC_flag_Tp = 1;
BC_flag_Tw = 1;

%flag_BC = 2; % 1: in vivo; 2: WinCF tube

SO_BC1 = 1; % Dirichlet BC of oxygen at x = L

theta_f_0 = 0.1; 
theta_p_0 = 0.2;

F_IC0 = 0;
P_IC0 = 0; 
I_IC0 = 0;  
SA_IC0 = 1;
SG_IC0 = 1;
SN_IC0 = 1;

Tf_IC0 = 2e-2; % total amount of antibiotic only kills anarobes added to the tube
Tp_IC0 = 2e-2; % total amount of antibiotic only kills Pseudomonas added to the tube
Tw_IC0 = 2e-2; % total amount of wide spectrum antibiotic added to the tube

Ind_Tp = 0;  % Indicator for antibiotic only kills Pseudomona
Ind_Tf = 0;  % Indicator for antibiotic only kills anarobes
Ind_Tw = 1;  % Indicator for wide spectrum antibiotic

%%%%%%%%%%%%%  domain information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 0.8;

Nx = 100;

dx = L/Nx;

t_Final = 50; % 100; % total simulation simulation time

dt = 0.01; % time step size

NT = floor(t_Final/dt);

Nstride = 100; % 10; %  Number of saved steps

FS_interval = floor(NT/Nstride); % File saving intervall



%N_Treat = 8;
% Treat_ID = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
% Treat_Names = ["NT", "Tp", "Tf", "Tw", "Tpf","Tpw","Tfw","Tpfw"];

N_Treat = 2;
Treat_ID = [0 0 0; 0 1 0];
Treat_Names = ["NT", "Tf"];

N_ratio = 1;
%Ratio_pf_ini = [0.2; 0.8];
Ratio_pf_ini = [0.7920 0.8397; 0.2052 0.1485];


%N_ratio = 24;
% Ratio_pf_ini = [0.0756  0.7920 0.1351 0.4841 0.8397 0.4549 0.1212 0.5900 0.7383 0.6596 0.5433 0.4201 0.7758 0.7007 0.2245 0.4564 0.8050 0.3693 0.9300 0.2317 0.6219 0.5562 0.7293 0.3323;
%                 0.9198  0.2052 0.8483 0.5114 0.1485 0.5288 0.8107 0.3846 0.2499 0.1025 0.4480 0.5666 0.2113 0.2852 0.5840 0.5378 0.1873 0.6089 0.0674 0.7605 0.2192 0.1264 0.2246 0.6673];
%             

Ratio_pf_ini = Ratio_pf_ini'; % make its dimension 24x2 instead of 2x24

% Sample_ID = [10 12 19 21 26 27];
% Sample_ID = [Sample_ID 31:48];

Sample_ID = [12 26];  % these two have biomass increase after antimicrobial treatments

Save_Inter_flag = 0; % 1: save intermediate values, 0: only save initial and final values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter to adjust pH inhibition on Pseudomonas, smaller lambda means
% stronger inhibition
lambda = 0.05; % 0.05: strong inhibition; 0.1: weak inhibition











