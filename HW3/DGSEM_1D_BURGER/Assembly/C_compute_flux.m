function  [flux] = C_compute_flux(femregion,u,Dati);
%%  [flux] = C_compute_flux(femregion,u,Dati);
%==========================================================================
% compute numerical flux according to the dg discretization with P0 el.
%==========================================================================
%    called in C_main1D_burger.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%          u           : (real)    solution vector
%
%    OUTPUT:
%          flux(u)     : godunov flux


addpath FESpace
addpath Assembly

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

% flux vector inizialization
flux = zeros(1,ndof);


for ie = 1 : ne
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
    [BJ, pphys_1D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_1D, femregion.h);
    
    
    x = pphys_1D;
    c2 = eval(Dati.c2);
    if (ie > 1)
        if(Dati.Burger == 0)
            flux(iglo) = flux(iglo) + c2(1)*u(iglo-1);
            flux(iglo) = flux(iglo) - c2(1)*u(iglo);
        else
            flux(iglo) = flux(iglo) + 0.5*c2(1)*u(iglo-1).^2;
            flux(iglo) = flux(iglo)- 0.5*c2(1)*u(iglo).^2;
        end
    else
        x = Dati.domain(1);
        t = Dati.t;
        phi = eval(Dati.exact_sol);
        if(Dati.Burger == 0)
            flux(iglo) = flux(iglo) + c2(1)*phi;
            flux(iglo) = flux(iglo) - c2(1)*u(iglo);
        else
            flux(iglo) = flux(iglo) + 0.5*c2(1)*phi^2;
            flux(iglo) = flux(iglo) - 0.5*c2(1)*u(iglo).^2;
        end
    end
    
end










