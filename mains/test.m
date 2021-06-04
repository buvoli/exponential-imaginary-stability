% -- amp functions --------------------------------------------------------
amp_sdc = @(q) @(z1, z2) abs(rESDC(z1, z2, struct('tau', lobpts(q, [0 1]), 'm', 2 * q - 3)));
amp_pbm = @(q) @(z1, z2) ampETDPBM(z1, z2, struct('z', [-1; legpts(q-1)], 'b', -1 * ones(q, 1), 'mS', 1, 'kappa', 1, 'alpha', 2));
amp_rk  = @(z1, z2) abs(rERK4(z1,z2));
%amp_rk  = @(z1, z2) abs(rIMRK4(z1,z2));
amps    = {amp_sdc(4), amp_pbm(5), amp_rk};

% -- angles ---------------------------------------------------------------
rhos = [0, pi/512, pi/256, pi/128 pi/56];
rho2angle = @(rho) 90 - atan(abs(real(exp(1i * (pi / 2 + rho))) / imag(exp(1i * (pi / 2 + rho))))) * 180 / pi;

% -- zs -------------------------------------------------------------------
Np = 200;
z1_real = linspace(0, 120, Np);
z2_real = linspace(-5, 5, Np);

%% == Generate Plots ======================================================

na = length(amps);
nr = length(rhos);

figure(1); clf;
figure_options = struct('FigureIndex', 1, 'ColoredLogAmp', true, 'ClearFigure', false, 'LogAmpCaxis', [-5,0]);

for i = 1 : na
    for j = 1 : nr
        subplot(nr, na, na * (j-1) + i);        
        z1_angle = exp(1i * (pi / 2 + rhos(j)));
        z2_angle = exp(1i * (pi / 2 - rhos(j)));        
        TDStabilityRegionPlot(amps{i}, z1_real, z2_real, z1_angle, z2_angle, figure_options)
    end
end

