% include directories
addpath(genpath('../../stepper/common'));
addpath(genpath('../../stepper/rk'));

% experiment settings
equation = 'zds';
epsilons = [0 0.5 1 2 4 6 8 10];
deltas = linspace(0, 8/100, 9);
integrator = @ERK;
integrator_params = struct('coeffGenerator', @ERK4, 'max_ts_to_store', 2);

%% == load equation =======================================================
wd = pwd();
cd(fullfile('../../stepper/equations/', equation))
init
cd(wd);

%% == compute reference using IMEXRK4 =====================================
fprintf('Computing reference... ');
pars.epsilon = 0;
[~, ys] = IMRK(LF(pars), NF, tspan, y0, Nt_reference, struct('coeffGenerator', @IMRK4, 'parameters', pars));
y_reference = ys(:,end);
fprintf('done.\n');

%% == run experiment with reference integrator ============================
num_Nts   = length(Nts);
imrk_error = zeros(num_Nts, 1);

pars.epsilon = 0; pars.delta = 0;
fprintf(' Running Reference Integrator = %i\n',  pars.epsilon);
for j = 1 : num_Nts
    fprintf('     Nt = %i/%i...', j , num_Nts);
    [~, ys] = IMRK(LF(pars), NF, tspan, y0, Nts(j), struct('coeffGenerator', @IMRK4, 'parameters', pars));
    imrk_error(j) = error_norm(ys(:,end), y_reference);
    fprintf('done.\n');
end

%% == run experiment with exponential integrator ==========================
% num_epsilons = length(epsilons);
% num_Nts = length(Nts);
% errors = zeros(num_epsilons, num_Nts);
% legend_labels = cell(num_epsilons, 1);
% 
% parfor i = 1 : num_epsilons
%     epars = mergeStructs(struct('epsilon', epsilons(i)), pars, false);
%     fprintf(' Running Epsilon = %i\n',  epars.epsilon);
%     for j = 1 : num_Nts
%         fprintf('     Nt = %i/%i...', j , num_Nts);
%         
%         [ts, ys] = integrator(LF(epars), NF, tspan, y0, Nts(j), mergeStructs(struct('parameters', epars), integrator_params, false));
%         errors(i,j) = error_norm(ys(:,end), y_reference);
%         
%         fprintf('done.\n');
%     end
%     legend_labels{i} = ['\epsilon = ', num2str(epsilons(i))];
% end

num_deltas = length(deltas);
num_Nts = length(Nts);
errors = zeros(num_deltas, num_Nts);
legend_labels = cell(num_deltas, 1);

parfor i = 1 : num_deltas
    epars = mergeStructs(struct('delta', deltas(i)), pars, false);
    fprintf(' Running delta = %i\n',  epars.epsilon);
    for j = 1 : num_Nts
        fprintf('     Nt = %i/%i...', j , num_Nts);
        
        [ts, ys] = integrator(LF(epars), NF, tspan, y0, Nts(j), mergeStructs(struct('parameters', epars), integrator_params, false));
        errors(i,j) = error_norm(ys(:,end), y_reference);
        
        fprintf('done.\n');
    end
    legend_labels{i} = ['\delta = ', num2str(deltas(i))];
end

%% == plot ================================================================
hs = diff(tspan) ./ Nts;

loglog(hs, errors);
legend(legend_labels);
xlabel('h');
ylabel('error');
axis([hs(end) hs(1) 1e-16 1e2])

hold on;
loglog(hs, imrk_error, 'k--');
hold off;

