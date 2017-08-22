%Shape function for 2D finite element method
%Input: xi, eta: natural coordinates, NPE: # of nodes per element
%       ELXY: coordinates of local nodes, array of size NPEx2
%       flag: 1 = linear triangle,  2 = quadratic triangle
%             3 = linear rectangle, 4 = quadratic rectangle
%Output: SF: value of shape functions, size NPEx1
%        GDSF: Derivate w.r.t. global coordinates, size NPEx2
%        DETJ: determinant of Jacobian
function [SF, GDSF, DETJ] = Shape2D(xi, eta, ELXY, NPE, flag)

% code common to all types of elements
SF = zeros(1,NPE);  % value of shape functions at (xi,eta)
GDSF = zeros(2,NPE); % derivatives w.r.t. global cooridinates
DSF = zeros(2,NPE);  % deirvatives w.r.t. natural coordinates
GJ = zeros(2);

if flag == 1 % linear triangle element
    
    SF = [xi, eta, 1 - xi - eta]; 
    DSF= [1, 0, -1; 0, 1, -1];
    
    
elseif flag == 2 % quadratic triangle element
        
    SF = [xi*(2*xi-1), eta*(2*eta-1), (1-xi-eta)*(1-2*xi-2*eta), ...
          4*xi*eta, 4*eta*(1-xi-eta), 4*xi*(1-xi-eta)];
    
    DSF = [4*xi-1, 0, 4*xi+4*eta-3, 4*eta, -4*eta, 4*(1-2*xi-eta); ...
           0, 4*eta-1, 4*xi+4*eta-3, 4*xi, 4*(1-xi-2*eta), -4*xi];
elseif flag == 3 % linear rectangle element
    
    SF = [(1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4];
    
    DSF = [-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4; ...
           -(1-xi)/4,  -(1+xi)/4, (1+xi)/4,  (1-xi)/4];
       
elseif flag == 4 % quadratic rectangle element
    
    SF = [(xi^2-xi)*(eta^2-eta)/4, (xi^2+xi)*(eta^2-eta)/4, (xi^2+xi)*(eta^2+eta)/4, (xi^2-xi)*(eta^2+eta)/4, ...
          (1-xi^2)*(eta^2-eta)/2, (xi^2+xi)*(1-eta^2)/2, (1-xi^2)*(eta^2+eta)/2, (xi^2-xi)*(1-eta^2)/2, (1-xi^2)*(1-eta^2)];
    
    DSF = [(2*xi-1)*(eta^2-eta)/4, (2*xi+1)*(eta^2-eta)/4, (2*xi+1)*(eta^2+eta)/4, (2*xi-1)*(eta^2+eta)/4, ...
           -xi*(eta^2-eta), (2*xi+1)*(1-eta^2)/2, -xi*(eta^2+eta), (2*xi-1)*(1-eta^2)/2, -2*xi*(1-eta^2);
           (xi^2-xi)*(2*eta-1)/4, (xi^2+xi)*(2*eta-1)/4, (xi^2+xi)*(2*eta+1)/4, (xi^2 - xi)*(2*eta+1)/4, ...
           (1-xi^2)*(2*eta-1)/2, -eta*(xi^2+xi), (1-xi^2)*(2*eta+1)/2, -eta*(xi^2-xi), -2*eta*(1-xi^2)];
end

    %Jacobian
    GJ = DSF*ELXY;
    DETJ = det(GJ);
    GDSF = inv(GJ)*DSF;