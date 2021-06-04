addpath(genpath('../../stability'));
addpath(genpath('../../stepper/common'));

% -- amp functions --------------------------------------------------------
amp_esdc = @(q) @(z1, z2) abs(rESDC(z1, z2, struct('tau', lobpts(q, [0 1]), 'm', 2 * q - 2)));
amp_epbm = @(q) @(z1, z2) rEPBM(z1, z2, struct('z', [-1; legpts(q-1)], 'b', -1 * ones(q, 1), 'mS', 1, 'kappa', 1, 'alpha', 1));
amp_erk  = @(z1, z2) abs(rERK4(z1,z2));
amp_imrk = @(z1, z2) abs(rIMRK4(z1,z2));

savename   = @(name, rho_index) fullfile('figures','rho', [name, '-rho-', num2str(rho_index)]);
savename3D = @(name, rho_index) fullfile('figures','vanilla', [name, '-3d']);
savenameV  = @(name) fullfile('figures','vanilla', [name, '-vanilla']);

Np = 501;
methods  = struct('amp', {amp_erk, amp_esdc(4), amp_epbm(5), amp_imrk}, 'name', {'erk', 'esdc', 'epbm', 'imrk'}, 'z2r', {linspace(-8,8,Np), linspace(-12,12,Np), linspace(-2,2, Np),  linspace(-6,6,Np)} );

amps     = {methods.amp};
m_names  = {methods.name};
z2_reals = {methods.z2r};

% -- angles ---------------------------------------------------------------
rhos = [0 pi/2048 pi/512 pi/128 pi/32];% [0, pi/256, pi/128 pi/64 pi/32];
rho2angle = @(rho) 90 - atan(abs(real(exp(1i * (pi / 2 + rho))) / imag(exp(1i * (pi / 2 + rho))))) * 180 / pi;

% -- zs -------------------------------------------------------------------
Np = 501;
z1_real = linspace(0, 60, Np);

%% == Generate Vanilla Stability Plots ====================================
na = length(amps);
fh = figure(100); clf;
figure_options = struct('FigureIndex', 100, 'ThreeDimensionalAmp', false, 'ClearFigure', true, 'FontSize', 15, 'LabelAxis', true, 'DrawAxis', false, 'ContourLevels', [0 1 1.01], 'ColorMap', [.4 .4 .4; .8 .8 .8]);

z1_angle = 1i;
z2_angle = 1i;

for i = 1 : na        
    TDStabilityRegionPlot(amps{i}, z1_real, z2_reals{i}, z1_angle, z2_angle, figure_options)
    hold on; plot([min(z1_real) max(z1_real)],[0 0], 'Color', [.3 .3 .3]); hold off; % fill in z=2 (otherwise it looks patchy)
    exportFigure(fh, struct('SavePath', savenameV(m_names{i}), 'Format', 'pdf', 'PaperPosition',   [0 0 8 8]))
end

%% == Generate Vanilla 3D Plots ===========================================
na = length(amps);
fh = figure(1); clf;
figure_options = struct('FigureIndex', 1, 'ThreeDimensionalAmp', true, 'ClearFigure', true, 'FontSize', 15, 'AmpCaxis', [.99, 1.01], 'ZLim', [.99, 1.01], 'LogAmp', false, 'ShowColorbar', false, 'DrawAxis', false, 'View', [53 45], 'ShowUnstable', true);

for i = 1 : na
    z1_angle = exp(1i * (pi / 2));
    z2_angle = exp(1i * (pi / 2)); 
    rho = 0;
    
    TDRStabilityRegionPlot(amps{i}, z1_real, z2_reals{i}, z1_angle, z2_angle, rho, figure_options)
    zlabel({'', '','$|R(iz_1,iz_2)|$'});        % Prevent export from clipping zlabel 
    xticks(z1_real([1 ceil(end/2) end]))        % limit x numticks to three
    yticks(z2_reals{i}([1 ceil(end/2) end]))    % limit y numticks to three
    exportFigure(fh, struct('SavePath', savename3D(m_names{i}, rho), 'Format', 'pdf', 'PaperPosition',   [0 0 8 8]));
    
end

%% == Generate 2D Rho Plots ===============================================
na = length(amps);
nr = length(rhos);

fh = figure(1); clf;
figure_options = struct('FigureIndex', 1, 'ThreeDimensionalAmp', true, 'ClearFigure', true, 'AmpCaxis', [-5,0], 'FontSize', 10, 'ShowColorbar', false, 'LabelAxis', false, 'DrawAxis', false);

z1_angle = 1i;
z2_angle = 1i;

for i = 1 : na
    for j = 1 : nr
        TDRStabilityRegionPlot(amps{i}, z1_real, z2_reals{i}, z1_angle, z2_angle, rhos(j), figure_options)
        exportFigure(fh, struct('SavePath', savename(m_names{i}, j), 'Format', 'pdf', 'PaperPosition',   [0 0 5 5]))
    end
end

%% --> save legend
fh = figure(1); clf;
colorbar;
caxis(figure_options.AmpCaxis);
set(gca, 'FontName', 'Minion Pro', 'FontSize', 10);
exportFigure(fh, struct('SavePath', fullfile('figures', 'rho', 'amp-colorbar'), 'Format', 'pdf', 'PaperPosition',   [0 0 5 5]))

%% == Generate 2D Rho Plots (Magnified -- show stability region separation) ========================================

methods  = struct('amp', {amp_erk, amp_esdc(4), amp_epbm(5), amp_imrk}, 'name', {'erk', 'esdc', 'epbm', 'imrk'}, 'z2r', {linspace(-3,3,Np), linspace(-3,3,Np), linspace(-3,3, Np), linspace(-1,1,Np)} );

amps     = {methods.amp};
m_names  = {methods.name};
z2_reals = {methods.z2r};
savenameRZ = @(name, rho_index) fullfile('figures','rho', [name, '-zoom-rho-', num2str(rho_index)]);

% -- angles ---------------------------------------------------------------
rhos = [0, pi/2048 pi/1048 pi/512 pi/256, pi/128 pi/64 pi/32 pi/16 pi/8 pi/6 pi/4 pi/3];

% -- zs -------------------------------------------------------------------
Np = 500;
z1_real = linspace(0, 6, Np);
z2_real = linspace(-1, 1, Np);

na = length(amps);
nr = length(rhos);

fh = figure(1); clf;
figure_options = struct('FigureIndex', 1, 'ThreeDimensionalAmp', true, 'ClearFigure', true, 'AmpCaxis', [-5,0], 'FontSize', 10, 'ShowColorbar', false, 'LabelAxis', false, 'DrawAxis', false);

z1_angle = 1i;
z2_angle = 1i;

for i = 1 : na
    for j = 1 : nr
        TDRStabilityRegionPlot(amps{i}, z1_real, z2_reals{i}, z1_angle, z2_angle, rhos(j), figure_options)
        exportFigure(fh, struct('SavePath', savenameRZ(m_names{i}, j), 'Format', 'pdf', 'PaperPosition',   [0 0 5 5]))
    end
end

% --> save legend
fh = figure(1); clf;
colorbar;
caxis(figure_options.AmpCaxis);
exportFigure(fh, struct('SavePath', fullfile('figures', 'exp', 'amp-colorbar'), 'Format', 'pdf', 'PaperPosition',   [0 0 5 5]))