%the script to solve a second order differential in one space dimension by
%finite element method


clear all;

t_start = cputime;

SetPara_WinCF_2D;

g0 = @(x) mu_f*(1-(2/pi)*atan(beta0*x));
g1 = @(x) K_f*(1 - 0.9*(2/pi)*atan(beta1*x));
g2 = @(x) K_p*(0.9*(0.5 - (1/pi)*atan(beta2*x)) + 0.1);

% the Gauss points and weights for quadrature on interval [-1 1]
global GAUSPT
GAUSPT = [0 0 0 0; -1/sqrt(3) 1/sqrt(3) 0 0; 0 0.7745966692 -0.7745966692 0; 0.33998104 -0.33998104 0.86113631 -0.86113631];

global GAUSWT
GAUSWT = [2 0 0 0; 1 1  0 0; 8/9 5/9 5/9 0; 0.65214515 0.65214515 0.34785485 0.34785485];

GAUSPT = GAUSPT';
GAUSWT = GAUSWT';

% the points and weights for quadrature on a triangle element (in
% barycentric coordinates
NTriQuad = 3;
TriQuadPT = zeros(NTriQuad,2);
TriQuadWT = zeros(NTriQuad,1);

if NTriQuad == 1
    TriQuadPT = [1/3 1/3];
    TriQuadWT = 1;
elseif NTriQuad == 3
    TriQuadPT = [1/2 0; 1/2 1/2; 0 1/2];
    TriQuadWT = [1/3; 1/3; 1/3];
elseif NTriQuad == 4
    TriQuadPT = [1/3 1/3; .6 .2; .2 .6; .2 .2];
    TriQuadWT = [-27/48; 25/48; 25/48; 25/48];
end

global NGP


NGP = 3;  % number of Gauss-points in numerical quadrature



flag = 3; % 1: linear triangle element, 3: linear rectangular element
imax = NY+1;

fid8 = fopen('./data/Total','w'); % file for total amount vs time



%%% ------------------------- Preprocessor -------------------------------


dx = (xR - xL)/NX; % mesh size in x-direction
dy = (yT - yB)/NY; % mesh size in y-direction

%preprocessor, generate mesh and connectivity matrix
[NNM, NEM, NPE, GLXY, NOD, NB_ELE_NO_F, NB_ELE_NO_SO, ...
    NB_ELE_F,  NB_ELE_SO, ...
    NB_EDG_F, NB_EDG_SO,  ...
    EB_NODE_NO_F,  EB_NODE_NO_SO,  ...
    EB_NODE_F,  EB_NODE_SO] = Mesh2D_WinCF(xL, xR, yB, yT, NX, NY, flag);

% initialize F, P, SO, SN concentration

[theta_f_n, theta_p_n, Fn, Pn, In, SOn, SNn, SAn, SGn, TBn, TTn] = Initialize_WinCF2D(NNM);

% dummy variables for the RHS function, only useful if Ind_B == 1 or Ind_T == 1

TBnph = zeros(size(Fn));
TBnp1 = TBnph;
TTnph = TBnph;
TTnp1 = TBnph;

Save_FPIONAGfp(Fn, Pn, In, SOn, SNn, SAn, SGn, TBn, TTn, theta_f_n, theta_p_n, imax, 0);
%calculate the elment equation and assemble the global matrix


ELXY = zeros(NPE,2);

Fn_EL = zeros(NPE,1);   % value of F  at time step n at nodes of the current element
Pn_EL = zeros(NPE,1);   % value of P  at time step n at nodes of the current element
In_EL = zeros(NPE,1);   % value of I  at time step n at nodes of the current element
SOn_EL = zeros(NPE,1);  % value of SO at time step n at nodes of the current element
SNn_EL = zeros(NPE,1);  % value of SN at time step n at nodes of the current element
SAn_EL = zeros(NPE,1);  % value of SA at time step n at nodes of the current element
SGn_EL = zeros(NPE,1);  % value of SG at time step n at nodes of the current element

if Ind_B == 1
    TBn_EL = zeros(NPE,1);  % value of TB at time step n at nodes of the current element
end
if Ind_T == 1
    TTn_EL = zeros(NPE,1);  % value of TT at time step n at nodes of the current element
end

theta_f_EL = zeros(NPE,1);   % value of fermenter cell density at nodes of current element
theta_p_EL = zeros(NPE,1);   % value of P. aeruginosa cell density at nodes of current element


%initialize, AF, AP, and AI, AB are independent of time, only initialize once
AF = spalloc(NNM,NNM,5*NNM);
AP = spalloc(NNM,NNM,5*NNM);
AI = spalloc(NNM,NNM,5*NNM);

if Ind_B == 1
    AB = spalloc(NNM,NNM,5*NNM);
end

flat_T = 1;
file_NO = 1;

for iT = 1:NT
    
     if iT*dt < 0.5*t_Final  % apply bicarbonate the first half
         Ind_B = 1;
         Ind_T = 0;
     else                    % apply tobramycin the second half
         Ind_B = 0;
         Ind_T = 1;
     end
         
     %%%%%%%%%%%%% Start Debugging %%%%%%%%%%%%%%%%%%
%       fprintf(1,'\n Ind = %d \n',Ind);
%       fprintf(1,'\n max(max(Fn - Pn)) = %6.5f, t_l = %6.5f, t_c = %6.5f \n', max(max(Fn - Pn)), t_l, t_c);
      %%%%%%%%%%%%% End Debugging %%%%%%%%%%%%%%%%%%%%
     
      % RHS for F
      FF = zeros(NNM,1); % the RHS vector
      
      % RHS for P
      FP = zeros(NNM,1); % the RHS vector
      
      % RHS for I
      FI = zeros(NNM,1); % the RHS vector
      
      % RHS for TB
      if Ind_B == 1
          FB = zeros(NNM,1); % the RHS vector
      end
      
      % matrix and RHS for SO
      AO = spalloc(NNM,NNM,5*NNM);
      FO = zeros(NNM,1); % the RHS vector
      
      % matrix and RHS for SA
      AN = spalloc(NNM,NNM,5*NNM);
      FN = zeros(NNM,1); % the RHS vector
      
      % matrix and RHS for SA
      AA = spalloc(NNM,NNM,5*NNM);
      FA = zeros(NNM,1); % the RHS vector
      
      % matrix and RHS for SG
      AG = spalloc(NNM,NNM,5*NNM);
      FG = zeros(NNM,1); % the RHS vector
      
      % matrix and RHS for TT
      if Ind_T == 1
          AT = spalloc(NNM,NNM,5*NNM);
          FT = zeros(NNM,1); % the RHS vector
      end
    
    % rectangle linear element
        
        for NE = 1:NEM %loop over all the elements
            
            for i = 1:NPE %get global coordinates of local nodes of element NE
                
                ELXY(i,1) = GLXY(NOD(NE,i),1);  % x-coordinate
                ELXY(i,2) = GLXY(NOD(NE,i),2);  % y-coordinate
                
                Fn_EL(i) = Fn(NOD(NE,i)); %
                Pn_EL(i) = Pn(NOD(NE,i)); %
                In_EL(i) = In(NOD(NE,i)); %
                SOn_EL(i) = SOn(NOD(NE,i)); %
                SNn_EL(i) = SNn(NOD(NE,i)); %
                SAn_EL(i) = SAn(NOD(NE,i)); %
                SGn_EL(i) = SGn(NOD(NE,i)); %
                
                if Ind_B == 1
                    TBn_EL(i) = TBn(NOD(NE,i)); %
                end
                
                if Ind_T == 1
                    TTn_EL(i) = TTn(NOD(NE,i)); %
                end
                
                
                theta_f_EL(i) = theta_f_n(NOD(NE,i));
                theta_p_EL(i) = theta_p_n(NOD(NE,i));
                
            end
            
            for NI = 1:NGP   %loop over Gauss point in x-direction
                for NJ = 1:NGP %loop over Gauss point in y-direction
                    
                    xi  = GAUSPT(NI,NGP);
                    eta = GAUSPT(NJ,NGP);
                    
                    [SF, GDSF, DETJ] = Shape2D(xi, eta, ELXY, NPE, flag);
                    
                    CONST = DETJ*GAUSWT(NI,NGP)*GAUSWT(NJ,NGP);
                    
                    xy = SF*ELXY;
                    
                    x = xy(1);
                    y = xy(2);
                    
                    Fn_Quad = SF*Fn_EL;  % value of F at time step n at the current quadrature point
                    Pn_Quad = SF*Pn_EL;
                    In_Quad = SF*In_EL;
                    SOn_Quad = SF*SOn_EL;
                    SNn_Quad = SF*SNn_EL;
                    SAn_Quad = SF*SAn_EL;
                    SGn_Quad = SF*SGn_EL;
                    
                    if Ind_B == 1
                        TBn_Quad = SF*TBn_EL; %
                    end
                    
                    if Ind_T == 1
                        TTn_Quad = SF*TTn_EL; %
                    end
                
                    theta_f_Quad = SF*theta_f_EL;
                    theta_p_Quad = SF*theta_p_EL;
                    
                    
                    A00_F = 1/dt;
                    A00_P = 1/dt;
                    A00_I = 1/dt;
                    
                    if Ind_B ==1
                        A00_TB = 1/dt;
                    end
                    
                    pHA = Fn_Quad - Pn_Quad + C_0;
                    
                    if Ind_B == 1
                        pHA = pHA - TBn_Quad;
                    end
                    
                    tmp_f = g0(SOn_Quad)*(1/(K_G + SGn_Quad))*theta_f_Quad*(1 - theta_f_Quad/g1(In_Quad));
                    tmp_pA = (mu_pa*SOn_Quad/(K_O + SOn_Quad) + exp(-d_p*SOn_Quad)*mu_pn*SNn_Quad/(K_N + SNn_Quad))*(1/(K_A + SAn_Quad))*theta_p_Quad*(1 - theta_p_Quad/g2(pHA));
                    tmp_pO = mu_pa*(1/(K_O + SOn_Quad))*(SAn_Quad/(K_A + SAn_Quad))*theta_p_Quad*(1 - theta_p_Quad/g2(pHA));
                    tmp_pN = exp(-d_p*SOn_Quad)*mu_pn*(1/(K_N + SNn_Quad))*(SAn_Quad/(K_A + SAn_Quad))*theta_p_Quad*(1 - theta_p_Quad/g2(pHA));
                    
                    A00_SO = 1/dt + tmp_pO/Y_pO + (d_O/Y_pO)*theta_p_Quad/(K_O + SOn_Quad);
                    A00_SN = 1/dt + tmp_pN/Y_pN;
                    A00_SA = 1/dt + tmp_pA/Y_pA;
                    A00_SG = 1/dt + tmp_f/Y_fG;
                    
                    
                    if Ind_T == 1
                        A00_TT = 1/dt + d_T*theta_p_Quad/Y_pT;
                    end
                    
                    
                    FXY_F =  Fn_Quad/dt + SGn_Quad*tmp_f/Y_fF;%  the right hand side for F
                    FXY_P =  Pn_Quad/dt + SOn_Quad*tmp_pO/Y_pPo + SNn_Quad*tmp_pN/Y_pPn; %    SAn_Quad*tmp_pA/Y_pP; %  the right hand side for P
                    FXY_I =  In_Quad/dt + SOn_Quad*tmp_pO/Y_pIo + SNn_Quad*tmp_pN/Y_pIn;%            SAn_Quad*tmp_pA/Y_pI; %  the right hand side for I
                    FXY_O =  SOn_Quad/dt; %  the right hand side for SO
                    FXY_N =  SNn_Quad/dt;  %  the right hand side for SN
                    FXY_A =  SAn_Quad/dt; %  the right hand side for SA
                    FXY_G =  SGn_Quad/dt; %  the right hand side for SG
                    
                    if Ind_B == 1
                        FXY_B = TBn_Quad/dt;
                    end
                    
                    if Ind_T == 1
                        FXY_T = TTn_Quad/dt;
                    end
                    %  FXY_N =  SOn_Quad/dt;
                    
                    for j = 1:NPE
                        
                        JJ = NOD(NE,j);
                        
                        FF(JJ) = FF(JJ) + CONST*SF(j)*FXY_F;
                        FP(JJ) = FP(JJ) + CONST*SF(j)*FXY_P;
                        FI(JJ) = FI(JJ) + CONST*SF(j)*FXY_I;
                        FO(JJ) = FO(JJ) + CONST*SF(j)*FXY_O;
                        FN(JJ) = FN(JJ) + CONST*SF(j)*FXY_N;
                        FA(JJ) = FA(JJ) + CONST*SF(j)*FXY_A;
                        FG(JJ) = FG(JJ) + CONST*SF(j)*FXY_G;
                        
                        
                        if Ind_B == 1
                            FB(JJ) = FB(JJ) + CONST*SF(j)*FXY_B;
                        end
                        
                        if Ind_T == 1
                            FT(JJ) = FT(JJ) + CONST*SF(j)*FXY_T;
                        end
                    
                        for i = 1:NPE
                            
                            II = NOD(NE,i);
                            
                            if flat_T == 1   % only calculate the matrix during the first time step, since it is independent of time
                                AF(II,JJ) = AF(II,JJ) + CONST*(D_F*GDSF(1,i)*GDSF(1,j) + D_F*GDSF(2,i)*GDSF(2,j) + A00_F*SF(i)*SF(j));
                                AP(II,JJ) = AP(II,JJ) + CONST*(D_P*GDSF(1,i)*GDSF(1,j) + D_P*GDSF(2,i)*GDSF(2,j) + A00_P*SF(i)*SF(j));
                                AI(II,JJ) = AI(II,JJ) + CONST*(D_I*GDSF(1,i)*GDSF(1,j) + D_I*GDSF(2,i)*GDSF(2,j) + A00_I*SF(i)*SF(j));
                                
                                
                                if Ind_B == 1
                                    AB(II,JJ) = AB(II,JJ) + CONST*(D_B*GDSF(1,i)*GDSF(1,j) + D_B*GDSF(2,i)*GDSF(2,j) + A00_TB*SF(i)*SF(j));
                                end
                                %flat_T = 2;
                            end
                            
                            AO(II,JJ) = AO(II,JJ) + CONST*(D_O*GDSF(1,i)*GDSF(1,j) + D_O*GDSF(2,i)*GDSF(2,j) + A00_SO*SF(i)*SF(j));
                            AN(II,JJ) = AN(II,JJ) + CONST*(D_N*GDSF(1,i)*GDSF(1,j) + D_N*GDSF(2,i)*GDSF(2,j) + A00_SN*SF(i)*SF(j));
                            AA(II,JJ) = AA(II,JJ) + CONST*(D_A*GDSF(1,i)*GDSF(1,j) + D_A*GDSF(2,i)*GDSF(2,j) + A00_SA*SF(i)*SF(j));
                            AG(II,JJ) = AG(II,JJ) + CONST*(D_G*GDSF(1,i)*GDSF(1,j) + D_G*GDSF(2,i)*GDSF(2,j) + A00_SG*SF(i)*SF(j));
                            
                            
                            if Ind_T == 1
                                AT(II,JJ) = AT(II,JJ) + CONST*(D_T*GDSF(1,i)*GDSF(1,j) + D_T*GDSF(2,i)*GDSF(2,j) + A00_TT*SF(i)*SF(j));
                            end
                            
                        end
                    end                   
                end              
            end          
        end      

    
    %Impose the boundary conditions, build condensed matrix equation,
    %here we have 4 chemical species, identified by the last argument
    
    
    % F, flag_spe = 1
    [AC_F, FC_F] = Boundary2D_WinCF(AF, FF, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 1);
    
    % P, flag_spe = 2
    [AC_P, FC_P] = Boundary2D_WinCF(AP, FP, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 2);
    
    % I, flag_spe = 3
    [AC_I, FC_I] = Boundary2D_WinCF(AI, FI, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 3);
    
    % SO, flag_spe = 3
    [AC_O, FC_O] = Boundary2D_WinCF(AO, FO, GLXY, NOD, NB_ELE_NO_SO, NB_ELE_SO, NB_EDG_SO, EB_NODE_NO_SO, EB_NODE_SO, flag, 4);
   % [AC_O, FC_O] = Boundary2D_WinCF(AO, FO, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 4);
   
    % SN, flag_spe = 7
    [AC_N, FC_N] = Boundary2D_WinCF(AN, FN, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 7);
    
    % SA, flag_spe = 5
    [AC_A, FC_A] = Boundary2D_WinCF(AA, FA, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 5);
    
    % SG, flag_spe = 6
    [AC_G, FC_G] = Boundary2D_WinCF(AG, FG, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 6);
    
    if Ind_B == 1
          [AC_B, FC_B] = Boundary2D_WinCF(AB, FB, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 1);
    end
    
    if Ind_T == 1
          [AC_T, FC_T] = Boundary2D_WinCF(AT, FT, GLXY, NOD, NB_ELE_NO_F, NB_ELE_F, NB_EDG_F, EB_NODE_NO_F, EB_NODE_F, flag, 1);
    end
    

    % Solve the condensed matrix equation

%     if(iT <= 3)
    
    Fnp1  = AC_F\FC_F;
    Pnp1  = AC_P\FC_P;
    Inp1  = AC_I\FC_I;   
    SOnp1 = AC_O\FC_O;
    SNnp1 = AC_N\FC_N;
    SAnp1 = AC_A\FC_A;
    SGnp1 = AC_G\FC_G;   
    
    if Ind_B == 1
         TBnp1 = AC_B\FC_B;
    end
    
    if Ind_T == 1
         TTnp1 = AC_T\FC_T;
    end
    
    
%     else
%         
%     Fnp1 = pcg(AC_F,FC_F,1e-5,200);
%     Pnp1 = pcg(AC_P,FC_P,1e-5,200);
%     SOnp1 = pcg(AC_O,FC_O,1e-5,200);
%     SNnp1 = pcg(AC_N,FC_N,1e-5,200);
% 
%     end
    
    Fnph = .5*(Fn + Fnp1);
    Pnph = .5*(Pn + Pnp1);
    Inph = .5*(In + Inp1);
    SOnph = .5*(SOn + SOnp1);
    SNnph = .5*(SNn + SNnp1);
    SAnph = .5*(SAn + SAnp1);
    SGnph = .5*(SGn + SGnp1);
    
    if Ind_B == 1
         TBnph = .5*(TBn + TBnp1);
    end
    
    if Ind_T == 1
         TTnph = .5*(TTn + TTnp1);
    end
    
    % use Runge Kutta method to solve the ODEs for theta_f and theta_p
    
    U = [theta_f_n, theta_p_n];
    
    Y1 = U;
    
    %Y = CF_FP_RHS_WinCF2D(U,SO,SA,SG,F,P,I)
    
    fY1 = CF_FP_RHS_WinCF2D(Y1,SOn,SNn,SAn,SGn,Fn,Pn,In,TBn,TTn);
    Y2 = U + .5*dt*fY1;
    
    fY2 = CF_FP_RHS_WinCF2D(Y2,SOnph,SNnph,SAnph,SGnph,Fnph,Pnph,Inph,TBnph,TTnph);
    Y3 = U + .5*dt*fY2;
    
    fY3 = CF_FP_RHS_WinCF2D(Y3,SOnph,SNnph,SAnph,SGnph,Fnph,Pnph,Inph,TBnph,TTnph);
    Y4 = U +  dt*fY3;
    
    fY4 = CF_FP_RHS_WinCF2D(Y4,SOnp1,SNnp1,SAnp1,SGnp1,Fnp1,Pnp1,Inp1,TBnp1,TTnp1);
    %update U
    U = U + (dt/6)*(fY1+ 2.0*fY2 + 2.0*fY3 + fY4);
    
    theta_f_n = U(:,1);
    theta_p_n = U(:,2);
    
    % update the chemicals
    Fn = Fnp1;
    Pn = Pnp1;
    In = Inp1;
    SOn = SOnp1;
    SNn = SNnp1;
    SAn = SAnp1;
    SGn = SGnp1;
      
    if Ind_B == 1
        TBn = TBnp1;
    end
    
    if Ind_T == 1
        TTn = TTnp1;
    end
    
    if  mod(iT,FS_interval) == 0
        
        %save data as 2D array
        %Save_FPONfp(Fn, Pn, SOn, SNn, theta_f_n, theta_p_n, imax, file_NO);
         
        Save_FPIONAGfp(Fn, Pn, In, SOn, SNn, SAn, SGn, TBn, TTn, theta_f_n, theta_p_n, imax, file_NO)
                
        fprintf(fid8, '%6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e', sum(Fn), sum(Pn), sum(In), sum(SOn), sum(SNn), sum(SAn), sum(SGn), sum(TBn), sum(TTn), sum(theta_f_n), sum(theta_p_n));
        fprintf(fid8, '\n');
        
%         %%%%%%%%%%%%% Start Debugging %%%%%%%%%%%%%%%%%%
%         fprintf(1,'\n Ind = %d \n',Ind);
%         fprintf(1,'\n max(max(Fn - Pn)) = %6.5f, t_l = %6.5f, t_c = %6.5f \n', max(max(Fn - Pn)), t_l, t_c);
%         %%%%%%%%%%%%% End Debugging %%%%%%%%%%%%%%%%%%%%
 
        
        file_NO = file_NO +1;

        fprintf(1,'\n Current time = %f \n', iT*dt);
        
        t_now = cputime;
        
        fprintf(1,'\n Time used  %8.6f seconds \n',t_now - t_start);
    end
    
    
    %hold off
    
    flat_T = 2;
    
    clear FF FP FI AO FO AN FN AA FA AG FG;
    
end


fclose(fid8);

% clf;
% 
[X, Y] = meshgrid(xL:dx:xR,yB:dy:yT);

t_end = cputime;

fprintf(1,'\n Total time used: %8.6f seconds \n',t_end - t_start);





