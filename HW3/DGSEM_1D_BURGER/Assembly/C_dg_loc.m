function [DG_loc] = C_dg_loc(dphiq,nln,c2)
%% [DG_loc] = C_dg_loc(dphiq,nln,c2)
%==========================================================================
% Build the local stiffness matrix for the term grad(u)grad(v)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%        dphiq         : basis functions
%          nln         : (integer) number of local unknowns
%          c2          : (array real) evaluation of wavespeed
%
%    OUTPUT:
%          K_loc       :  (array real) Local stiffness matrix


DG_loc = zeros(nln,nln);

%% General implementation -- to be used with general finite element spaces
for i=1:nln
    for j=1:nln
        
        for k = 1
            DG_loc(i,j) = DG_loc(i,j) +  c2(k) * dphiq(1,k,j) * dphiq(1,k,i);
        end
        
    end
end



                                              
                                              

