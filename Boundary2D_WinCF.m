%function to impose boundary conditions (Essential and Natural) for 2D finite element code

% input flag_spe : 1 -> F, 2 -> P, 3 -> SO, 4 -> SN 
function [AC, FC] = Boundary2D_WinCF(A, F, GLXY, NOD, NB_ELE_NO, NB_ELE, NB_EDG, EB_NODE_NO, EB_NODE, flag, flag_spe)


AC = A;
FC = F;

global GAUSPT     
global GAUSWT
global NGP

% edge-node connectivity 
Lin_TRI_EN = [1 2; 2 3; 3 1];  % linear triangle
Quad_TRI_EN = [1 4 2; 2 5 3; 3 6 1]; % quadratic rectangle
Lin_REC_EN = [1 2; 2 3; 3 4; 4 1]; % linear rectangle
Quad_REC_EN = [1 5 2; 2 6 3; 3 7 4; 4 8 1]; % quadratic rectangle

% Natural boundary condition

for IBEL = 1:NB_ELE_NO % loop over natural boundary element
    
    GBE = NB_ELE(IBEL,1); % global element number
    
    for IBED = 1:NB_ELE(IBEL,2) % loop over boundary edges 
        
        LBED = NB_EDG(IBEL,IBED); % local edge # of the boundary edge
        
        if flag == 1
            LN1 = Lin_TRI_EN(LBED,1);
            LN2 = Lin_TRI_EN(LBED,2);

        elseif flag == 2
            LN1 = Quad_TRI_EN(LBED,1);
            LN2 = Quad_TRI_EN(LBED,3);

        elseif flag == 3
            LN1 = Lin_REC_EN(LBED,1);
            LN2 = Lin_REC_EN(LBED,2);

        elseif flag == 4
            LN1 = Quad_REC_EN(LBED,1);
            LN2 = Quad_REC_EN(LBED,3);
        end
        
        x1 = GLXY(NOD(GBE,LN1),1);
        y1 = GLXY(NOD(GBE,LN1),2);
        
        x2 = GLXY(NOD(GBE,LN2),1);
        y2 = GLXY(NOD(GBE,LN2),2);
        
        h = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        for NI = 1:NGP
            
            xi = GAUSPT(NI,NGP);
            CONST = GAUSWT(NI,NGP)*h/2;  %Jacobian is h/2
            
            if flag == 1 || flag == 3 % linear element
                SF = SF1D(xi,1);
                s = SF*[0, h]';
            end
            
            x = (1 - s/h)*x1 + (s/h)*x2;
            y = (1 - s/h)*y1 + (s/h)*y2;
            
            if flag == 1
                for j = 1:2  %loop over node on this edge
                    J = NOD(GBE,Lin_TRI_EN(LBED,j)); % global node number of this local node
                    FC(J) = FC(J) + SF(j)*g_N(x,y,flag_spe)*CONST;
                end  
            
            elseif flag == 3
                 for j = 1:2  %loop over node on this edge
                    J = NOD(GBE,Lin_REC_EN(LBED,j)); % global node number of this local node
                    FC(J) = FC(J) + SF(j)*g_N(x,y,flag_spe)*CONST;
                end
            end   
        end
    end
end
            
% Essential boundary condition

for i = 1:EB_NODE_NO % loop over Essential boundary node
    I = EB_NODE(i);
    x = GLXY(I,1);
    y = GLXY(I,2);
    FC = FC - AC(:,I)*g_E(x,y,flag_spe);
    FC(I) = g_E(x,y,flag_spe);
    AC(:,I) = 0;
    AC(I,:) = 0;
    AC(I,I) = 1;
end


% essential boundary condition
function z = g_E(x,y,flag_s)


global  SO_BC0 

if flag_s == 4 % SO
    
    z = SO_BC0;
    
else
    
    z = 0;
    
end

return

% natural boundary condition
function z = g_N(x,y,flag_s)

  z = 0; % homogeneous NBC for all chemical species

return


% function return 1-D shape function at natural coordinate xi, order 1 or 2
function z = SF1D(xi,order)

if order == 1
    z = [(1-xi)/2, (1+xi)/2];
elseif order == 2
    z = [xi*(xi-1)/2, -(xi-1)*(xi+1), xi*(xi+1)/2];
end

return



