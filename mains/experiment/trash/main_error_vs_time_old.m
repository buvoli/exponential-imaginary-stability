% include directories
addpath(genpath('../../stepper/common'));
addpath(genpath('../../stepper/rk'));
addpath(genpath('../../stepper/poly'));
addpath(genpath('../../stepper/sdc'));

%% == load equation =======================================================
equation = 'kdv';
num_times = 100;

wd = pwd();
cd(fullfile('../../stepper/equations/', equation))
init
cd(wd);

%% experiment settings
integrator = 'esdc'; % erk, esdc, epbm
partitioning = 'rotation'; % rotation, rotation-uxx, translation, hyperdiffusion

par_to_label_value = @(v) num2str(v);
par_to_label_name  = @(v) ['\', v];
switch( partitioning )
    case 'rotation'
        par_name = 'rho';
        par_values = [0 pi/256 pi/128 pi/64 pi/32];  
        par_to_label_value = @rhoLabel; 
    case 'rotation-uxx'
        par_name = 'rhouxx';
        par_values = [0 pi/256 pi/128 pi/64 pi/32];
        par_to_label_value = @rhoLabeAlex Garrido Outonl;
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
    
    otherwise
        error('invalid partitioning');
end

switch ( integrator )
    case 'erk'
        integrator = @ERK;
        integrator_options = struct('coeffGenerator', @ERK4);
        line_marker = 's-';
    case 'esdc'
        integrator = @ESDC;
        integrator_options = struct('tau', lobpts(4, [0 1]), 'm', 6); 
        line_marker = 'd-';
    case 'epbm'
        integrator = @EPBM_PMFCmS;
        integrator_options = struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 0, 'alpha', 1);
        %integrator_options = struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 1, 'alpha', 1);
        %integrator_options = struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 2, 'alpha', 2);
        %integrator_options = struct('z', lobpts(4), 'b', -1 * ones(4), 'mS', [], 'kappa', 1, 'alpha', 2);
        line_marker = '*-';
    otherwise
        error('invalid integrator');
end

%% == compute reference using RK4 =========================================
fprintf('Computing reference... ');
reference_ppars = mergeStructs(struct('epsilon', 0), pars, false);
reference_options = struct('parameters', reference_ppars, 'max_ts_to_store', num_times);
L_ref = LF(reference_ppars);
RHS = @(t,y,pars) L_ref .* y + NF(t,y,pars);
[~, ~, ts_ref, ys_ref] = rk4(RHS, tspan, y0(:), Nt_reference, reference_options);
% reference_options = struct('coeffGenerator', @IMRK4, 'parameters', reference_ppars);
% [~, ys] = IMRK(LF(reference_ppars), NF, tspan, y0, Nt_reference, reference_options);
% y_reference = ys(:,end);
fprintf('done.\n');

%% == prepare run parameters ==============================================
num_pars = length(par_values);
runs     = cell(num_pars, 1);
for i = 1 : num_pars
    prb_pars = mergeStructs(struct(par_name, par_values(i)), pars, false); % problem parameters
    runs{i} = struct( ...
        'integrator', integrator, ...
        'options', mergeStructs(struct('parameters', prb_pars, 'max_ts_to_store', num_times), integrator_options), ...
        'legend', ['$', par_to_label_name(par_name), ' = ', par_to_label_value(par_values(i)), '$']);
end
epbm
% -- add IMRK4 integrator ---
imrk4_options = struct('coeffGenerator', @IMRK4, 'parameters', reference_ppars, 'max_ts_to_store', num_times);
runs{end+1} = struct('integrator', @IMRK, 'options', imrk4_options, 'legend', 'IMRK4');

%% == run integrators =====================================================
num_ts   = size(ys_ref,2);
num_runs = length(runs);
errors = zeros(num_ts, num_runs);
parfor i = 1 : num_runsepbm
    prob_pars = runs{i}.options.parameters;
    fprintf(' Running %s = %i\n',  par_name, prob_pars.(par_name));
    [ts, ys] = runs{i}.integrator(LF(prob_pars), NF, tspan, y0, Nt_run, runs{i}.options);
    errors(:, i) = error_norm(ys, ys_ref);
    fprintf('done.\n');
end

%% == produce plots =======================================================
f = figure();

line_markers = {'-', '-<', '->', '-^', '-v'};

MS = 10 - 5 * (0:num_pars - 1) / (num_pars - 1); MS = [MS(end) MS(1:end-1)];
LW = 3.5 - 1 * (0:num_pars - 1) / (num_pars - 1); LW = [LW(end) LW(1:end-1)];

for i = 1 : length(par_values)
    h = semilogy(ts_ref, errors(:, i), line_markers{i}, 'LineWidth', LW(i), 'MarkerSize', MS(i)); hold on;
    set(h, 'MarkerFaceColor', min(1, get(h, 'Color') * 1))
    set(h, 'MarkerEdgeColor', min(1, get(h, 'Color') * 1))
end
hold off;


xlabel('time $t$', 'Interpreter', 'Latex');
ylabel('Relative Error','Interpreter', 'Latex');
axis([ts_ref(1) ts_ref(end) 1e-8 1e2])

% place rho = 0 line on top.
chH = get(gca,'Children');
set(gca,'Children',[chH(end); chH(1:end-1)]);

hold on;
semilogy(ts_ref, errors(:, end), 'k--', 'LineWidth', 2);
hold off;

%xticks([1e-3, 1e-2, 1e-1])
yticks(10.^(-8:2:2))
%legend(legend_labels, 'Location', 'NorthWest', 'Interpreter', 'latex'); legend boxoff;

f_info = functions(integrator);
exportFigure(f, struct('Format', 'pdf', 'SavePath', fullfile('output', [equation, '-time-', num2str(Nx), '-', f_info.function, '-', partitioning]), 'PaperPosition', [0 0 8 7]));


%% Legend Figure
% f = figure();
% 
% legend_labels = cellfun(@(run) run.legend, runs, 'UniformOutput', false);
% for i = 1 : length(par_values)
%     h = plot(0, 0, line_markers{i}, 'LineWidth', LW(i), 'MarkerSize', MS(i));
%     set(h, 'MarkerFaceColor', min(1, get(h, 'Color') * 1))
%     set(h, 'MarkerEdgeColor', min(1, get(h, 'Color') * 1))
%     hold on;
% end
% legend(legend_labels, 'Location', 'northoutside', 'Interpreter', 'latex', 'Orientation', 'horizontal'); legend boxoff;
% set(gca, 'FontSize', 12);
% exportFigure(f, struct('Format', 'pdf', 'SavePath', fullfile('output', ['legend-', partitioning]), 'PaperPosition', [0 0 20 8]));

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
    l = [num2str(omega), '\hspace{1em}'];
end