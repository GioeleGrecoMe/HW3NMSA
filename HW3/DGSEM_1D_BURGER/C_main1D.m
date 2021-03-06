function [errors,solutions,femregion,Dati] = C_main1D(TestName,h)
%==========================================================================
% Solution of the Wave Equation with linear finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          TestName    : (string)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Dati        : (struct)  see C_dati.m
%
% Usage:
%    [errors,solutions,femregion,Dati] = C_main1D('Test1',3)



addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing
addpath SemLib


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati(TestName);
Ltot=Dati.domain;
Li=Ltot(1);
Lf=Ltot(2);
nRef=round(log2((Lf-Li)/h));
Dati.nRefinement = nRef;

%==========================================================================
% MESH GENERATION
%==========================================================================

% [Region] = C_create_mesh(Dati);
[Region] = C_create_mesh_sem(Dati);


%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================
x = Dati.domain(1);
if  eval(Dati.c2) > 0
    Dati.InflowPoint(1) = 1;
else
    Dati.InflowPoint(1) = 0;
end

[M_nbc,A_nbc] = C_matrix1D(Dati,femregion);


%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================
Dati.t = 0;
t = 0;

%==========================================================================
% BUILD INITIAL CONDITIONS
%==========================================================================
x = femregion.coord;
u0 = eval(Dati.u0);
[u0] = C_snapshot_1D(femregion,u0,Dati);

%==========================================================================
% STARTING THE TIME-LOOP
%==========================================================================

fprintf('============================================================\n')
fprintf('Starting time-loop ... \n');
fprintf('============================================================\n')
time_surf = [0: Dati.dt : Dati.T];

% FOR PLOTTING u(x,t) -- store soultions ---
u_surf = zeros(length(time_surf),size(u0,1));
u_surf(1,:) = u0;
k_surf = 2;


for t = 0 : Dati.dt : Dati.T - Dati.dt
    
    fprintf('time = %5.3e \n',t);
    
    Dati.t = t;
    % Mu' = flux(u);
    % M(u1 - u0) = dt*flux(u0)
    [flux] = C_compute_flux(femregion,u0,Dati);
    b_nbc = M_nbc*u0 + Dati.dt*flux';
    u1 = M_nbc\b_nbc;

    % M(u2 - (3/4u0+1/4u1)) = 1/4dt*flux(u1) 
    Dati.t = t + Dati.dt;
    [flux] = C_compute_flux(femregion,u1,Dati);
    b_nbc = 0.25*M_nbc*u1 +0.75*M_nbc*u0 + 0.25*Dati.dt*flux';
    u2 = M_nbc\b_nbc;
    
    % M(u3 - (1/3*u0 + 2/3*u2)) = 2/3*dt*flux(u2)
    Dati.t = t + 0.5*Dati.dt;
    [flux] = C_compute_flux(femregion,u2,Dati);
    b_nbc = M_nbc*(1/3*u0 + 2/3*u2) + 2/3*Dati.dt*flux';
    u3 =  M_nbc\b_nbc;
    [u3] = C_snapshot_1D(femregion,u3,Dati);
    
    u_surf(k_surf,:) = u3;
    k_surf = k_surf + 1;
    
    u0 = u3;
    
end

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================
[solutions] = C_postprocessing(Dati,femregion,u1);

figure(100);
time = 0 : Dati.dt : Dati.T;
surf(x,time, u_surf,'EdgeColor','None'); hold on;
gf = gca;
xlim([min(x)  max(x)]); ylim([0 Dati.T]);
gf.XTick = [min(x)  max(x)];
gf.XTickLabel = {num2str(Dati.domain(1)), num2str(Dati.domain(2))};
gf.YTick = [0 Dati.T/2 Dati.T];
gf.YTickLabel = {num2str(0),num2str(Dati.T*0.5) ,num2str(Dati.T)};
view(2); xlabel('space-axis'); ylabel('time-axis'); title('u_h(x,t)');



%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
% if (Dati.plot_errors)
%     [errors] = C_compute_errors(Dati,femregion,solutions);
% end



