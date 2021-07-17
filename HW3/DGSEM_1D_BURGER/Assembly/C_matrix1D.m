function [M,A]=C_matrix1D(Dati,femregion)
%% [M,A] = C_matrix1D(Dati,femregion)
%==========================================================================
% Assembly of the mass matrix M and the stiffness matrix A
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          M           : (sparse(ndof,ndof) real) mass matrix
%          A           : (sparse(ndof,ndof) real) stiffnes matrix


addpath FESpace
addpath Assembly

fprintf('============================================================\n')
fprintf('Assembling matrices M and A ... \n');
fprintf('============================================================\n')


% connectivity infos
ndof         = femregion.ndof; % degrees of freedom
nln          = femregion.nln;  % local degrees of freedom
ne           = femregion.ne;   % number of elements
connectivity = femregion.connectivity; % connectivity matrix


if nln > 1
    % quadrature nodes and weights for integrals
    [nodes_1D,w_1D] = xwlgl(nln);
    
    % evaluation of shape bases and their derivative
    [dphiq,Grad] = basis_and_der_at_lgl(nodes_1D,nln);
else
    nodes_1D = 0;
    w_1D =  2;
    dphiq(:,:,1) = 1;
    Grad(:,:,1) = [0 0]';
    
end

% Assembly begin ...
M = sparse(ndof,ndof);  % Global Mass matrix
A = sparse(ndof,ndof);  % Global Stiffness matrix

for ie = 1 : ne
    
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
    
    [BJ, pphys_1D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_1D, femregion.h);
    
    % BJ        = Jacobian of the elemental map
    % pphys_2D = vertex coordinates in the physical domain
    
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    
    % Local stiffness matrix
    % ATT: Now c^2 has to be evaluated inside the integral!!!
    x = pphys_1D; c2 = eval(Dati.c2);
    
    [Adv_loc] = C_adv_loc(Grad,dphiq,w_1D,nln,BJ,c2);
    A(iglo,iglo) = A(iglo,iglo) + Adv_loc;
    
    
   
    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    
    % Local mass matrix
    [M_loc] = C_mass_loc(dphiq,w_1D,nln,BJ);
    
    % Assembly phase for mass matrix
    M(iglo,iglo) = M(iglo,iglo) + M_loc;
    
    
    %=============================================================%
    % DG MATRIX
    %=============================================================%
    [DG_loc] = C_dg_loc(dphiq,nln,c2);
    A(iglo,iglo) = A(iglo,iglo) + DG_loc;
    
    %     if(ie < ne & && Dati.InflowPoint(1) ~= 1)
    %
    %         %for the neighbouring element ie+1
    %         iglo_ne = connectivity(1:nln,ie+1);
    %         x = pphys_1D_ne; c2_ne = eval(Dati.c2);
    %
    %         [DG_loc_ne] = C_dg_loc_ne(dphiq,c2_ne,nln);
    %         A(iglo,iglo_ne) = A(iglo,iglo_ne) + DG_loc_ne;
    %     end
    
    if(ie > 1 && Dati.InflowPoint(1) == 1)
        
        %for the neighbouring element ie-1
        iglo_ne = connectivity(1:nln,ie-1);
        c2_ne = eval(Dati.c2);
        
        [DG_loc_ne_1] = C_dg_loc_ne_1(dphiq,c2_ne,nln);
        A(iglo,iglo_ne) = A(iglo,iglo_ne) + DG_loc_ne_1;
        
    end
    
    
    
    
    
    
    
end
