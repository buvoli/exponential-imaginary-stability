% include directories
addpath(genpath('../../stepper/common'));
addpath(genpath('../../stepper/rk'));
addpath(genpath('../../stepper/poly'));
addpath(genpath('../../stepper/sdc'));

% experiment settings
equation = 'zds';

%% == load equation =======================================================
wd = pwd();
cd(fullfile('../../stepper/equations/', equation))
init
cd(wd);

%% == compute reference using RK4 =====================================
fprintf('Computing reference... ');
reference_ppars = mergeStructs(struct('rho', 0, 'epsilon', 0), pars, false);
reference_options = struct('parameters', reference_ppars, 'coeffGenerator', @RK4C);
L_ref = LF(reference_ppars);
RHS = @(t,y,pars) L_ref .* y + NF(t,y,pars);
[~, y_reference] = RK(RHS, tspan, y0(:), Nt_reference, reference_options);
fprintf('done.\n');

%% == run integrators =====================================================

runs = {
    struct('integrator', @ERK,          'legend', 'ERK4',  'options', struct('coeffGenerator', @ERK4, 'parameters', pars), 'linestyle', 's-')
    struct('integrator', @ESDC,         'legend', 'ESDC6', 'options', struct('tau', lobpts(4, [0 1]), 'm', 5, 'parameters', pars), 'linestyle', 'd-')
    struct('integrator', @EPBM_PMFCmS,  'legend', 'EPBM5', 'options', struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 1, 'parameters', pars, 'alpha', 1), 'linestyle', 'o-')
    struct('integrator', @IMRK,         'legend', 'IMRK4', 'options', struct('coeffGenerator', @IMRK4, 'parameters', pars), 'linestyle', 'k*--')
    struct('integrator', @(L,N, tspan, y0, Nt, options) rk4(RHS, tspan, y0(:), Nt, options), 'legend', 'RK4', 'options', struct('parameters', pars), 'linestyle', 'p-')
};

num_Nts  = length(Nts);
num_runs = length(runs);
errors = zeros(num_Nts, num_runs);
parfor i = 1 : num_runs
    prob_pars = runs{i}.options.parameters;
    fprintf(' Running Integrator %s\n', runs{i}.legend);
    for j = 1 : num_Nts
        
        
        prob_pars.stepsize = tspan(end) / Nts(j);
        fprintf('     Nt = %i/%i...', j , num_Nts);
        
        [ts, ys] = runs{i}.integrator(LF(prob_pars), NF, tspan, y0, Nts(j), runs{i}.options);
        errors(j, i) = error_norm(ys(:,end), y_reference);
        
        fprintf('done.\n');
    end
end

save(['data-vanilla-', equation, '.mat'], 'runs', 'tspan', 'Nts', 'errors', 'equation')

%% == produce plot =======================================================
f = figure();

hs = diff(tspan) ./ Nts;
legend_labels = cellfun(@(run) run.legend, runs, 'UniformOutput', false);

for i = 1 : length(runs)
    h = loglog(hs, errors(:, i), runs{i}.linestyle, 'LineWidth', 3, 'MarkerSize', 10); hold on;
    set(h, 'MarkerFaceColor', min(1, get(h, 'Color') * 1.5))
end

xlabel('Stepsize $h$', 'Interpreter', 'Latex');
ylabel('Relative Error','Interpreter', 'Latex');
axis([hs(end) hs(1) 1e-10 1e2])
hold off;
set(gca, 'FontSize', 14, 'FontName', 'Minion Pro');

legend(legend_labels, 'Location', 'NorthWest', 'Interpreter', 'latex'); legend boxoff;
exportFigure(f, struct('SavePath', ['output/vanilla_integrators_', equation], 'Format', 'pdf', 'PaperPosition', [0 0 15 12], 'Renderer', 'opengl'));