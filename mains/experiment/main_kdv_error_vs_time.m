% include directories
addpath(genpath('../../stepper/common'));
addpath(genpath('../../stepper/rk'));
addpath(genpath('../../stepper/poly'));
addpath(genpath('../../stepper/sdc'));

%% == load equation =======================================================
equation = 'kdv';
num_times = 25; % note must be divisor of NT_reference and NT_run (or ts will not line up)

wd = pwd();
cd(fullfile('../../stepper/equations/', equation))
init
cd(wd);

if(mod(Nt_reference, num_times) ~= 0 || mod(Nt_run, num_times) ~= 0)
    error('num_times must be divisor of NT_reference and NT_run (or ts will not line up)'); 
end

%% experiment settings
partitioning = 'translation'; % rotation, rotation-uxx, translation, hyperdiffusion

parameters_erk  = struct('blank', 0, 'rho', pi/64, 'rhouxx', pi/3, 'epsilon', 16);
parameters_esdc = struct('blank', 0, 'rho', pi/64, 'rhouxx', pi/3, 'epsilon', 16);
parameters_epbm = struct('blank', 0, 'rho', pi/64, 'rhouxx', pi/3, 'epsilon', 16);
parameters      = fieldnames(parameters_erk);
par_legend_label = struct('rho', @rhoLabel, 'rhouxx', @rhoUxxLabel, 'epsilon', @epsilonLabel, 'blank', @blankLabel);

colors = {[0, 0.4470, 0.7410], .2 * [255 255 255]/255, .5 * [255 255 255]/255, .8 * [255 255 255]/255};

% integrator options
erk_integrator_options  = struct('coeffGenerator', @ERK4);
esdc_integrator_options = struct('tau', lobpts(4, [0 1]), 'm', 6); 
epbm_integrator_options = struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 1, 'alpha', 1);

%% == compute reference using RK4 =========================================
fprintf('Computing reference... ');
reference_ppars = mergeStructs(struct('epsilon', 0, 'rho', 0, 'rhouxx', 0), pars, false);
reference_options = struct('parameters', reference_ppars, 'max_ts_to_store', num_times);
L_ref = LF(reference_ppars);
RHS = @(t,y,pars) L_ref .* y + NF(t,y,pars);
[~, ~, ts_ref, ys_ref] = RK(RHS, tspan, y0(:), Nt_reference, mergeStructs(reference_options, struct('coeffGenerator', @RK4C)));
fprintf('done.\n');

%% == prepare run parameters ==============================================
num_pars = length(parameters);
runs     = cell(3 * num_pars, 1);
for i = 1 : num_pars
    
    par_name = parameters{i};
    
    par_value = parameters_erk.(par_name);
    prb_pars = mergeStructs(struct(par_name, par_value), pars, false); % problem parameters
    runs{3 * (i - 1) + 1} = struct( ...
        'integrator', @ERK, ...
        'options', mergeStructs(struct('parameters', prb_pars, 'max_ts_to_store', num_times), erk_integrator_options), ...
        'linemarker', 's-', ...
        'legend', ['$', par_legend_label.(par_name)(par_value), '$']);
    
    par_value = parameters_esdc.(par_name);
    prb_pars = mergeStructs(struct(par_name, parameters_esdc.(par_name)), pars, false); % problem parameters
    runs{3 * (i - 1) + 2} = struct( ...
        'integrator', @ESDC, ...
        'options', mergeStructs(struct('parameters', prb_pars, 'max_ts_to_store', num_times), esdc_integrator_options), ...
        'linemarker', 'd-', ...    
        'legend', ['$', par_legend_label.(par_name)(par_value), '$']);

    par_value = parameters_epbm.(par_name);
    prb_pars = mergeStructs(struct(par_name, parameters_epbm.(par_name)), pars, false); % problem parameters
    runs{3 * (i - 1) + 3} = struct( ...
        'integrator', @EPBM_PMFCmS, ...
        'options', mergeStructs(struct('parameters', prb_pars, 'max_ts_to_store', num_times), epbm_integrator_options), ...
        'linemarker', '*-', ...
        'legend', ['$', par_legend_label.(par_name)(par_value), '$']);
    
end
% % -- add IMRK4 integrator ---
% imrk4_options = struct('coeffGenerator', @IMRK4, 'parameters', reference_ppars, 'max_ts_to_store', num_times);
% runs{end+1} = struct('integrator', @IMRK, 'options', imrk4_options, 'legend', 'IMRK4');

%% == run integrators =====================================================
num_ts   = size(ys_ref,2);
num_runs = length(runs);
errors = zeros(num_ts, num_runs);
parfor i = 1 : num_runs
    prob_pars = runs{i}.options.parameters;
    fprintf(' Running integrator %i/%i\n',  i, num_runs);
    [ts, ys] = runs{i}.integrator(LF(prob_pars), NF, tspan, y0, Nt_run, runs{i}.options);
    errors(:, i) = error_norm(ys, ys_ref);
    fprintf('done.\n');
end

save(['data-long-time-', equation, '-t', num2str(tspan(end)), '-Nx', num2str(Nx)], 'errors', 'runs', 'num_pars', 'ts_ref', 'colors', 'equation', 'Nx');

%% == produce plots =======================================================
f = figure();

line_markers = {'-', '-<', '->', '-^', '-v'};
line_style = {'-', '-', '-.', ':'};

MS = 12 - 5 * (0:num_pars - 1) / (num_pars - 1); MS = [MS(end) MS(1:end-1)];
LW = 8 - 6 * (0:num_pars - 1) / (num_pars - 1); LW = [LW(end) LW(1:end-1)];

MS = 11 + 0 * (0:num_pars - 1) / (num_pars - 1); MS = [MS(end) MS(1:end-1)];
LW = 2 + 2 * (0:num_pars - 1) / (num_pars - 1); LW = [LW(end) LW(1:end-1)];

for i = 1 : num_pars
        
    mi = i+1:num_pars:length(errors); % marker indices
    ind = 3 * (i - 1) + 1;    
    if(i == 1)
        mi = union(mi, find(~isnan(errors(:, ind)), 1, 'last')); % add last non nan
    end
    
    h = semilogy(ts_ref, errors(:, ind), ['s', line_style{i}], 'LineWidth', LW(i), 'MarkerSize', MS(i), 'Color', colors{i}, 'MarkerIndices', mi); hold on;
    set(h, 'MarkerFaceColor', min(1, get(h, 'Color') * 1))
    set(h, 'MarkerEdgeColor', min(1, get(h, 'Color') * .9))
    
    mi = i+1:num_pars:length(errors); % marker indices
    ind = 3 * (i - 1) + 2;    
    if(i == 1)
        mi = union(mi, find(~isnan(errors(:, ind)), 1, 'last')); % add last non nan
    end
    
    h = semilogy(ts_ref, errors(:, ind), ['d', line_style{i}], 'LineWidth', LW(i), 'MarkerSize', MS(i), 'Color', colors{i}, 'MarkerIndices', mi); hold on;
    set(h, 'MarkerFaceColor', min(1, get(h, 'Color') * 1))
    set(h, 'MarkerEdgeColor', min(1, get(h, 'Color') * .9))
    
    mi = i+1:num_pars:length(errors); % marker indices
    ind = 3 * (i - 1) + 3;
    if(i == 1)
        mi = union(mi, find(~isnan(errors(:, ind)), 1, 'last')); % add last non nan
    end
    
    h = semilogy(ts_ref, errors(:, ind), ['o', line_style{i}], 'LineWidth', LW(i), 'MarkerSize', 0.8 * MS(i), 'Color', colors{i}, 'MarkerIndices', mi); hold on;
    set(h, 'MarkerFaceColor', min(1, get(h, 'Color') * 1))
    set(h, 'MarkerEdgeColor', min(1, get(h, 'Color') * .9))
end
hold off;

xlabel('Time $t$', 'Interpreter', 'Latex');
ylabel('Relative Error','Interpreter', 'Latex');
axis([ts_ref(1) ts_ref(end) 1e-6 1e4])

% place rho = 0 line on top.
chH = get(gca,'Children');
set(gca,'Children',[chH([end, end-1, end-2]); chH(1:end-3)]);

yticks(10.^(-8:2:2))
set(gca, 'FontSize', 14, 'FontName', 'Minion Pro');

exportFigure(f, struct('Format', 'pdf', 'SavePath', fullfile('output', [equation, '-', num2str(Nx)], [equation, '-error-T', num2str(tspan(end))]), 'PaperPosition', [0 0 8*3 7*1.5]));


%% Legend Figure
f = figure();

legend_labels = cellfun(@(run) run.legend, runs, 'UniformOutput', false);
eqn_dir = [equation, '-', num2str(Nx)];

% -- repartition legend 
for i = 1 : num_pars
    h = plot(0, 0, [line_style{i}], 'LineWidth', 3, 'MarkerSize', MS(i), 'Color', colors{i}, 'MarkerFaceColor', min(1, colors{i}), 'MarkerEdgeColor', colors{i} * .9); hold on;
end
legend(legend_labels(1:3:end), 'Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex'); legend boxoff;
set(gca, 'FontSize', 12);
exportFigure(f, struct('Format', 'pdf', 'SavePath', fullfile('output', eqn_dir , 'ltc-legend-repartition'), 'PaperPosition', [0 0 6*num_pars 4]));

% -- method legend 
f = figure();
plot(0, 0, 's', 'LineWidth', 3, 'MarkerSize', MS(i), 'Color', 'k', 'MarkerFaceColor',[1 1 1], 'MarkerEdgeColor', [.3 .3 .3]); hold on;
plot(0, 0, 'd', 'LineWidth', 3, 'MarkerSize', MS(i), 'Color', 'k', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',  [.3 .3 .3]); hold on;
plot(0, 0, 'o', 'LineWidth', 3, 'MarkerSize', MS(i), 'Color', 'k', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',  [.3 .3 .3]); hold on;
legend({'ERK4', 'ESDC6', 'EPBM4'}, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex'); legend boxoff;
set(gca, 'FontSize', 12);
exportFigure(f, struct('Format', 'pdf', 'SavePath', fullfile('output', eqn_dir , 'ltc-legend-integrator'), 'PaperPosition', [0 0 6*num_pars 4]));

function l = rhoLabel(rho)
    if(rho == 0)
        l = '\epsilon = 0';
    else
        [~, b] = rat(rho / pi);
        l = ['\mathbf{D} = -\mathrm{diag}(|\mathbf{k}|^3),~  \rho = \frac{\pi}{', num2str(b),'}']; %\mathbf{D} = \text{diag}(i|\mathbf{k}|^3), 
    end
    l = [l, '\hspace{1em}'];
end

function l = rhoUxxLabel(rho)
    if(rho == 0)
        l = '\mathbf{D} = \rho = 0';
    else
        [a, b] = rat(rho / pi);
        if(a == 1)       
            l = ['\mathbf{D} = -\mathrm{diag}(\mathbf{k}^2),~ \rho = \frac{\pi}{', num2str(b),'}'];
        else
            l = ['\mathbf{D} = -\mathrm{diag}(\mathbf{k}^2),~ \rho = \frac{', num2str(a), '\pi}{', num2str(b),'}'];
        end   
    end
    l = [l, '\hspace{1em}'];
end

function l = epsilonLabel(epsilon)
    l = ['\mathbf{D} = -\mathbf{I},~ \epsilon = ', num2str(epsilon)];
    l = [l, '\hspace{1em}'];
end

function l = blankLabel(rho)
    l = rhoLabel(0);
end