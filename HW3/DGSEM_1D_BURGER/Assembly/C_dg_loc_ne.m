function [DG_loc] = C_dg_loc_ne(Grad,dphiq,BJ,c2,...
                                Grad_ne,BJ_ne,c2_ne,nln,teta, penalty)
%% [DG_loc] = C_dg_loc(Grad,w_1D,nln,BJ)
%==========================================================================
% Build the local stiffness matrix for the term grad(u)grad(v)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Grad        : (array real) evaluation of the gradient on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%          c2          : (array real) evaluation of wavespeed
%
%    OUTPUT:
%          K_loc       :  (array real) Local stiffness matrix


DG_loc = zeros(nln,nln);

%% General implementation -- to be used with general finite element spaces
for i=1:nln
    for j=1:nln
        
        Binv_ne = 1./BJ_ne;    % inverse
        Binv = 1./BJ;    % inverse
        
        DG_loc(i,j) = DG_loc(i,j) -  0.5 * c2_ne(1) * dphiq(1,nln,i) * (Grad_ne(1,:,j) * Binv_ne )' ...
                                  - teta*0.5 * c2(nln) * dphiq(1,1,j) * (Grad(nln,:,i) * Binv )' ...
                                  - penalty * dphiq(1,1,j) * dphiq(1,nln,i);
        
    end
end



                                              
                                              

