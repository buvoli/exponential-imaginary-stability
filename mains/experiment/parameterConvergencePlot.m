function parameterConvergencePlot(equation, integrator, partitioning)
%PARAMETERCONVERGENCEPLOT Summary of this function goes here
%   Detailed explanation goes here

settings = getSettings(equation, integrator, partitioning);
[errors, hs, settings] = getData(equation, partitioning, settings);
makeConvergePlot(equation, partitioning, errors, hs, settings);
makeLegend(partitioning, settings);

end

function settings = getSettings(equation, integrator, partitioning)

% Load Equation       
switch( equation )
    case {'zds', 'kdv'}
    otherwise
        error('invalid equation');
end

% Parameter Settings
par_to_label_name  = @(v) ['\', v];
switch( partitioning )
    case 'rotation'
        par_name = 'rho';
        par_values = [0 pi/2048 pi/512 pi/128 pi/32]; % [0 pi/256 pi/128 pi/64 pi/32];
        par_to_label_value = @rhoLabel;
    case 'rotation-uxx'
        par_name = 'rhouxx';
        par_values = [0 pi/256 pi/64 pi/16 pi/4]; %[0 pi/256 pi/128 pi/64 pi/32];
        par_to_label_value = @rhoLabel;
        par_to_label_name  = @(v) '\rho';
    case 'translation'
        par_name = 'epsilon';
        par_values = [0 1 2 4 8];
        par_to_label_value = @epsilonLabel;
    case {'hyperdiffusion-2', 'hyperdiffusion-4', 'hyperdiffusion-6', 'hyperdiffusion-8'}
    
        par_name   = 'hyperv_coeff';
        par_to_label_name  = @(v) '\omega';
        par_to_label_value = @omegaLabel;
        
        switch(partitioning)
            case 'hyperdiffusion-2'
                pars.hyperv_order = 2;
                par_values = [0 1e6 1e8 1e10 1e12];
            case 'hyperdiffusion-4'
                pars.hyperv_order = 4;
                par_values = [0 1e4 1e6 1e8 1e10];
            case 'hyperdiffusion-6'
                pars.hyperv_order = 6;
                par_values = [0 1e2 1e4 1e6 1e8];
            case 'hyperdiffusion-8'
                pars.hyperv_order = 8;
                par_values = [0 1 1e2 1e4 1e6];
        end
        
        if(integrator ~= "erk")
            warning('hyperviscosity parameters are only tuned for ERK4');
        end
    
    otherwise
        error('invalid partitioning');
end

switch ( integrator )
    case 'erk'
        integrator = @ERK;
        integrator_options = struct('coeffGenerator', @ERK4, 'max_ts_to_store', 2);
        integrator_linemarker = 's-';
    case 'esdc'
        integrator = @ESDC;
        integrator_options = struct('tau', lobpts(4, [0 1]), 'm', 6); 
        integrator_linemarker = 'd-';
    case 'epbm'
        integrator = @EPBM_PMFCmS;
        integrator_options = struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 1, 'alpha', 1);
        integrator_linemarker = 'o-';
    otherwise
        error('invalid integrator');
end

settings = struct();
settings.par_name = par_name;
settings.par_values = par_values;
settings.par_to_label_name = par_to_label_name;
settings.par_to_label_value = par_to_label_value;
settings.integrator = integrator;
settings.integrator_options = integrator_options;
settings.integrator_linemarker = integrator_linemarker;

end


function [errors, hs, settings] = getData(equation, partitioning, settings)

y0 = []; tspan = []; Nts = []; LF = []; NF = []; error_norm = []; Nx = [];  % blank variables (overriten by load; required for parfor)
%% load equation        
wd = pwd();
cd(fullfile('../../stepper/equations/', equation))
init
cd(wd);

switch(partitioning) 
    case 'hyperdiffusion-2'
        pars.hyperv_order = 2;
    case 'hyperdiffusion-4'
        pars.hyperv_order = 4;
    case 'hyperdiffusion-6'
        pars.hyperv_order = 6;
    case 'hyperdiffusion-8'
        pars.hyperv_order = 8;
end

%% compute reference using RK4
fprintf('Computing reference... ');
reference_ppars = mergeStructs(struct('epsilon', 0, 'rho', 0, 'rhouxx', 0), pars, false);
reference_options = struct('parameters', reference_ppars, 'coeffGenerator', @RK4C);
L_ref = LF(reference_ppars);
RHS = @(t,y,pars) L_ref .* y + NF(t,y,pars);
[~, y_reference] = RK(RHS, tspan, y0(:), Nt_reference, reference_options);
fprintf('done.\n');

%% prepare run parameters
num_pars = length(settings.par_values);
runs     = cell(num_pars, 1);
for i = 1 : num_pars
    prb_pars = mergeStructs(struct(settings.par_name, settings.par_values(i)), pars, false); % problem parameters
    runs{i} = struct( ...
        'integrator', settings.integrator, ...
        'options', mergeStructs(struct('parameters', prb_pars), settings.integrator_options), ...
        'legend', ['$', settings.par_to_label_name(settings.par_name), ' = ', settings.par_to_label_value(settings.par_values(i)), '$']);
end

% -- add IMRK4 integrator ---
imrk4_options = struct('coeffGenerator', @IMRK4, 'parameters', reference_ppars);
runs{end+1} = struct('integrator', @IMRK, 'options', imrk4_options, 'legend', 'IMRK4');

%% == run integrators =====================================================
num_Nts  = length(Nts);
num_runs = length(runs);
errors = zeros(num_Nts, num_runs);
parfor i = 1 : num_runs
    prob_pars = runs{i}.options.parameters;
    fprintf(' Running %s = %i\n',  settings.par_name, prob_pars.(settings.par_name));
    for j = 1 : num_Nts
        
        prob_pars.method_stepsize = tspan(end) / Nts(j);
        fprintf('     Nt = %i/%i...', j , num_Nts);
        
        [~, ys] = runs{i}.integrator(LF(prob_pars), NF, tspan, y0, Nts(j), runs{i}.options);
        errors(j, i) = error_norm(ys(:,end), y_reference);
        
        fprintf('done.\n');
    end
end

hs = diff(tspan) ./ Nts;
settings.legend_labels = cellfun(@(run) run.legend, runs, 'UniformOutput', false);
settings.Nx = Nx;

end

% == Start Plotting Functions =============================================

function makeConvergePlot(equation, partitioning, errors, hs, settings)

num_pars = length(settings.par_values);

lM = lineMarkers();
lW = lineWidths(num_pars);
mS = markerSizes(num_pars);

f = figure();
for i = 1 : num_pars
    h = loglog(hs, errors(:, i), lM{i}, 'LineWidth', lW(i), 'MarkerSize', mS(i)); hold on;
    color = get(h, 'Color');
    set(h, 'MarkerFaceColor', color)
    set(h, 'MarkerEdgeColor', color)
end
hold off;

xlabel('Stepsize $h$', 'Interpreter', 'Latex');
ylabel('Relative Error','Interpreter', 'Latex');
axis([hs(end) hs(1) 1e-10 1.2e2])

% place rho = 0 line on top.
chH = get(gca,'Children');
set(gca,'Children',[chH(end); chH(1:end-1)]);

hold on;
loglog(hs, errors(:, end), 'k--', 'LineWidth', 2);
hold off;

xticks([1e-3, 1e-2, 1e-1])
yticks(10.^(-10:2:2))

f_info = functions(settings.integrator);
exportFigure(f, struct('Format', 'pdf', 'SavePath', fullfile('output', [equation, '-', num2str(settings.Nx), '-', f_info.function, '-', partitioning]), 'PaperPosition', [0 0 8 7]));

end

function makeLegend(partitioning, settings)

num_pars = length(settings.par_values);

lM = lineMarkers();
lW = lineWidths(num_pars);
mS = markerSizes(num_pars);

f = figure();

for i = 1 : num_pars
    h = plot(0, 0, lM{i}, 'LineWidth', lW(i), 'MarkerSize', mS(i));
    color = get(h, 'Color');
    set(h, 'MarkerFaceColor', color)
    set(h, 'MarkerEdgeColor', color)
    hold on;
end

plot(0, 0, 'k--', 'LineWidth', 2);
hold off;

legend(settings.legend_labels, 'Location', 'northoutside', 'Interpreter', 'latex', 'Orientation', 'horizontal'); legend boxoff;
set(gca, 'FontSize', 12);
exportFigure(f, struct('Format', 'pdf', 'SavePath', fullfile('output', ['legend-', partitioning]), 'PaperPosition', [0 0 25 8]));

end

function line_markers = lineMarkers()
    line_markers = {'-', '-<', '->', '-^', '-v'};
end

function marker_size = markerSizes(num_pars)
    marker_size = 10 - 5 * (0:num_pars - 1) / (num_pars - 1); 
    marker_size = [marker_size(end) marker_size(1:end-1)];
end

function line_width = lineWidths(num_pars)
    line_width = 3.5 - 1 * (0:num_pars - 1) / (num_pars - 1);
    line_width = [line_width(end) line_width(1:end-1)];
end

% == Start Label Functions ================================================

function l = rhoLabel(rho)
    if(rho == 0)
        l = '0 \hspace{1em}';
    else
        [~, b] = rat(rho / pi);
        l = ['\pi / ', num2str(b), '\hspace{1em}'];
    end
end

function l = epsilonLabel(epsilon)
    l = [num2str(epsilon), '\hspace{1em}'];
end

function l = omegaLabel(omega)
    if(omega == 0)
        l = ['0\hspace{1em}'];
        return;
    end
    
    pow = floor( log(omega) / log(10) + 2 * eps );
    mnt = omega * 10^(-pow);
    if(mnt == 1)
        l = ['10^{', num2str(pow), '}\hspace{1em}'];
    else    
        l = [num2str(mnt), '\times 10^{', num2str(pow), '}\hspace{1em}'];
    end
end
