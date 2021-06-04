function cs = L(params)
%L Returns coefficients for linear ODEs
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
ks = [0:Nx/2 -Nx/2+1:-1].' * (2*pi/Lx);

delta_rot     = abs(ks).^3 / tan(pi/2 + params.rho);                        % D = diag(abs(ks.^3))
delta_rot_uxx = ks.^2 / tan(pi/2 + params.rhouxx);                          % D = diag(abs(ks.^2))
delta_trans   = -params.epsilon;                                            % D = -I
delta_total   = params.delta^2 * (delta_rot + delta_rot_uxx) + delta_trans;

cs = -params.delta^2 * (1i*ks).^3  + delta_total; %- params.epsilon_tr - params.epsilon_rt * abs(ks).^5;
if(params.antialias)
    cs = antialias(cs);
end
end

