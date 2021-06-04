function cs = L(params)
%L Returns coefficients for linear ODEs
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
ks = [0:Nx/2 -Nx/2+1:-1].' * (2*pi/Lx);

delta_rot     = epsilon(params.rho) * abs(ks).^3;                           % D = diag(abs(ks.^3))
delta_rot_uxx = epsilon(params.rhouxx) * ks.^2;                                % D = diag(abs(ks.^2))
delta_trans   = - params.epsilon;                                           % D = -I
delta_hypr    = - (params.method_stepsize)^params.method_order * params.hyperv_coeff * ks.^params.hyperv_order; % hyperviscosity
delta_total   = delta_rot + delta_rot_uxx + delta_trans + delta_hypr;

cs = (1i*ks).^3 + delta_total;
if(params.antialias)
    cs = antialias(cs);
end
end

function e = epsilon(rho)
    if(rho == 0)
        e = 0;
    else
        e = 1 / tan(pi/2 + rho);
    end
end