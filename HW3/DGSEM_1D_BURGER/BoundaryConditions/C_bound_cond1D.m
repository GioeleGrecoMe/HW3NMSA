function [A,b] = C_bound_cond1D(A,b,femregion,Dati)
%% [A,b] = C_bound_cond1D(A,b,femregion,Dati)
%==========================================================================
% Assign Dirchlet boundary conditions
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%          femregion   : (struct)  see C_create_femregion.m
%          Dati        : (struct)  see C_dati.m

%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffness matrix
%          b           : (sparse(ndof,1) real) rhs vector
%

%fprintf('============================================================\n')
%fprintf('Assign Dirichlet boundary conditions ... \n');
%fprintf('============================================================\n')


ndof = length(b);
u_g = sparse(ndof,1);

boundary_points = femregion.boundary_points;
x = femregion.dof(boundary_points,1);
t = Dati.t;
u_g(boundary_points) = eval(Dati.exact_sol); % Compute the lifting operator ug



% Reduce the system A in order to solve the pb with 
% homogeneous Dirichlet conditions 
for k = 1:length(boundary_points)
    if(boundary_points(k) == Dati.InflowPoint) 
       A(boundary_points(k),:) = 0;
       A(boundary_points(k),boundary_points(k)) = 1;
       b(boundary_points(k)) = u_g(boundary_points(k));
    end
end


