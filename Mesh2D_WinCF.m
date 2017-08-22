%mesh generation function of a rectangular domain for 2D finite element
%method
%Input: xL, xR (x-coordinates of left and right boundary)
%       yB, yT (y-coordinates of bottom and top boundary)
%       NX, NY (number of intervals in x and y direction)
%       flag = 1: linear triangle, 2: quadratic triangle,
%              3: linear rectangle, 4: quadratic rectangle
%       flag_bc = 1: Essential BC at left and top boundaries
%                 2: EBC at left boundary
%                 3: EBC at top boundary
%Output: NNM: # of nodes in the mesh, NEM: # of elements in the mesh
%        NPE: # of nodes per element
%        GLXY: coordinates of global nodes, array of size NNMx2
%        GLXY(I,1) = x-coord of I-th global node, GLXY(I,2) = y-coord
%        NOD: connectivity matrix, array of size NEMxNPE
%        NB_ELE_NO: number of elements with Natural BC
%        NB_ELE: element list with natural BC, array of size NB_ELE_NOx2
%        NB_ELE(I,1): element index of I-th Natural boundary element
%        NB_ELE(I,2): number of natural boundary condition edges with NBC in I-th boundary element
%        NB_EDG: local edge number of the boundary element, array of size
%                NB_ELE_NOx2, NB_EDG(I,j) = local edge number of j-th
%                boundary edge in I-the boundary element, 1 <= j <= NBELE(I,2)
%        EB_NODE_NO: number of nodes with Essential BCs
%        EB_NODE: global node number of EBC node, array of size
%                 EB_NODE_NOx1
% Note F, P, I, SA, SG have same BCs, No flux BC at all boundary
% SO have no flux at left, right, bottom, and Dirichlet at top

function [NNM, NEM, NPE, GLXY, NOD, NB_ELE_NO_F, NB_ELE_NO_SO, ...
    NB_ELE_F,  NB_ELE_SO,  ...
    NB_EDG_F,  NB_EDG_SO,  ...
    EB_NODE_NO_F,  EB_NODE_NO_SO,  ...
    EB_NODE_F,  EB_NODE_SO] = Mesh2D_WinCF(xL, xR, yB, yT, NX, NY, flag)

dx = (xR - xL)/NX;
dy = (yT - yB)/NY;

if flag == 1 % linear triangle
    
    % size of various structures
    NNM = (NY+1)*(NX+1);
    NEM = 2*NX*NY;
    NPE = 3;
    
    GLXY = zeros(NNM,2);
    NOD = zeros(NEM,NPE);
    
    % global coordinate of nodes
    for I = 1:NNM
        GLXY(I,1) = mod(I-1,NX+1)*dx + xL;
        GLXY(I,2) = floor((I-1)/(NX+1))*dy + yB;
    end
    
    % connectivity matrix
    for j = 1:NY
        for i = 1:NX
            
            ELN1 = (j-1)*2*NX + 2*i - 1; %lower element number
            ELN2 = ELN1 + 1;             %upper element number
            
            NOD(ELN1,1) = (j-1)*(NX+1) + i;
            NOD(ELN1,2) = NOD(ELN1,1) + 1;
            NOD(ELN1,3) = NOD(ELN1,2) + NX + 1;
            
            NOD(ELN2,1) = NOD(ELN1,1);
            NOD(ELN2,2) = NOD(ELN1,3);
            NOD(ELN2,3) = NOD(ELN2,2) - 1;
            
        end
    end
    
    % BC for F  P, I, SA, SG
    
    % no-flux at all boundaries
    
    NB_ELE_NO_F = 2*NX + 2*NY;
    NB_ELE_F = zeros(NB_ELE_NO_F,2);
    NB_EDG_F = zeros(NB_ELE_NO_F,1);
    
    for i = 1:NX 
        % bottom boundary
        NB_ELE_F(i,1) = 2*i-1; 
        NB_ELE_F(i,2) = 1;
        NB_EDG_F(i,1) = 1;
    
      
        % top boundary
        NB_ELE_F(NX+i,1) = 2*NX*(NY-1) + 2*i;    % global element #
        NB_ELE_F(NX+i,2) = 1;          % # of edges with natural BC in i-th boundary element
        NB_EDG_F(NX+i,1) = 2;
    end
    
    for j = 1:NY
        % left boundary
        NB_ELE_F(2*NX+j,1) = (j-1)*2*NX+2; 
        NB_ELE_F(2*NX+j,2) = 1;
        NB_EDG_F(2*NX+j,1) = 3;
    
        % right boundary
        NB_ELE_F(2*NX+NY+j,1) = j*2*NX-1; 
        NB_ELE_F(2*NX+NY+j,2) = 1;
        NB_EDG_F(2*NX+NY+j,1) = 2;        
    end
    
    % Essential at left, right, bottom
    
    EB_NODE_NO_F = 0;  % # of nodes with EBC
    EB_NODE_F = zeros(EB_NODE_NO_F,1);
    
    % BC for SO, no flux at left, right, bottom, and Dirichlet at top
    
    NB_ELE_NO_SO = NX + 2*NY;
    NB_ELE_SO = zeros(NB_ELE_NO_SO,2);
    NB_EDG_SO = zeros(NB_ELE_NO_SO,1);
    
     for i = 1:NX 
        % bottom boundary
        NB_ELE_SO(i,1) = 2*i-1; 
        NB_ELE_SO(i,2) = 1;
        NB_EDG_SO(i,1) = 1;      
     end
     
    for j = 1:NY
        % left boundary
        NB_ELE_SO(NX+j,1) = (j-1)*2*NX+2; 
        NB_ELE_SO(NX+j,2) = 1;
        NB_EDG_SO(NX+j,1) = 3;
    
        % right boundary
        NB_ELE_SO(NX+NY+j,1) = j*2*NX-1; 
        NB_ELE_SO(NX+NY+j,2) = 1;
        NB_EDG_SO(NX+NY+j,1) = 2;        
    end
    % Essential boundary condition at the entire boundary
    
    EB_NODE_NO_SO = NX+1;  % # of nodes with EBC
    EB_NODE_SO = zeros(EB_NODE_NO_SO,1);
    
    %top boundary
    for i = 1:NX+1
        EB_NODE_SO(i) = NY*(NX+1) + i;
    end
    
   elseif flag == 3 % linear rectangle element

    % size of various structures
    NNM = (NX+1)*(NY+1);
    NEM = NX*NY;
    NPE = 4;
        
    GLXY = zeros(NNM,2);
    NOD = zeros(NEM,NPE);
    
    % global coordinate of nodes
    for I = 1:NNM
        GLXY(I,1) = mod(I-1,NX+1)*dx + xL;
        GLXY(I,2) = floor((I-1)/(NX+1))*dy + yB;
    end

        % connectivity matrix
    for j = 1:NY
        for i = 1:NX
            
            ELN = (j-1)*NX + i; % element number
            
            NOD(ELN,1) = (j-1)*(NX+1) + i;
            NOD(ELN,2) = NOD(ELN,1) + 1;
            NOD(ELN,3) = NOD(ELN,2) + NX + 1;
            NOD(ELN,4) = NOD(ELN,3) - 1;
            
        end
    end

     % BC for F  P, I, SA, SG
    
    % no-flux at all boundaries
    
    NB_ELE_NO_F = 2*NX + 2*NY;
    NB_ELE_F = zeros(NB_ELE_NO_F,2);
    NB_EDG_F = zeros(NB_ELE_NO_F,1);
    
    for i = 1:NX 
        % bottom boundary
        NB_ELE_F(i,1) = i; 
        NB_ELE_F(i,2) = 1;
        NB_EDG_F(i,1) = 1;
    
      
        % top boundary
        NB_ELE_F(NX+i,1) = NX*(NY-1) + i;    % global element #
        NB_ELE_F(NX+i,2) = 1;          % # of edges with natural BC in i-th boundary element
        NB_EDG_F(NX+i,1) = 3;
    end
    
    for j = 1:NY
        % left boundary
        NB_ELE_F(2*NX+j,1) = (j-1)*NX+1; 
        NB_ELE_F(2*NX+j,2) = 1;
        NB_EDG_F(2*NX+j,1) = 4;
    
        % right boundary
        NB_ELE_F(2*NX+NY+j,1) = j*NX; 
        NB_ELE_F(2*NX+NY+j,2) = 1;
        NB_EDG_F(2*NX+NY+j,1) = 2;        
    end
    
    % Essential at left, right, bottom
    
    EB_NODE_NO_F = 0;  % # of nodes with EBC
    EB_NODE_F = zeros(EB_NODE_NO_F,1);
    
    % BC for SO, no flux at left, right, bottom, and Dirichlet at top
    
    NB_ELE_NO_SO = NX + 2*NY;
    NB_ELE_SO = zeros(NB_ELE_NO_SO,2);
    NB_EDG_SO = zeros(NB_ELE_NO_SO,1);
    
     for i = 1:NX 
        % bottom boundary
        NB_ELE_SO(i,1) = i; 
        NB_ELE_SO(i,2) = 1;
        NB_EDG_SO(i,1) = 1;      
     end
     
    for j = 1:NY
        % left boundary
        NB_ELE_SO(NX+j,1) = (j-1)*NX+1;
        NB_ELE_SO(NX+j,2) = 1;
        NB_EDG_SO(NX+j,1) = 4;
    
        % right boundary
        NB_ELE_SO(NX+NY+j,1) = j*NX; 
        NB_ELE_SO(NX+NY+j,2) = 1;
        NB_EDG_SO(NX+NY+j,1) = 2;        
    end
    % Essential boundary condition at the entire boundary
    
    EB_NODE_NO_SO = NX+1;  % # of nodes with EBC
    EB_NODE_SO = zeros(EB_NODE_NO_SO,1);
    
    %top boundary
    for i = 1:NX+1
        EB_NODE_SO(i) = NY*(NX+1) + i;
    end
    
else
    
    fprintf(1,'\n flag must be 1 or 3, only use linear triangular or rectangular element for now. \n');
    
end














