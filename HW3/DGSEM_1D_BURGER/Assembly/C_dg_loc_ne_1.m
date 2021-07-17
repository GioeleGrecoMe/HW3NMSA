function [DG_loc] = C_dg_loc_ne_1(dphiq,c2_ne,nln)
%% [DG_loc] = C_dg_loc_ne_1(dphiq,c2_ne,nln)
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
        
        DG_loc(i,j) = DG_loc(i,j) - c2_ne(nln) * dphiq(1,nln,j) * dphiq(1,1,i);
        
    end
end



                                              
                                              

