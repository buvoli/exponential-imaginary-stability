% =========================================================================
% Main file for ZDS Repartitioning and Hyperviscosity Experiments
% =========================================================================

% include directories
addpath(genpath('../../stepper/common'));
addpath(genpath('../../stepper/rk'));
addpath(genpath('../../stepper/poly'));
addpath(genpath('../../stepper/sdc'));

equation = 'zds';

%% Repartitioning Experiments
integrators = {'erk', 'esdc', 'epbm'};
partitionings = {'rotation', 'rotation-uxx', 'translation'};

for i = 1 : length(integrators)
    for j = 1 : length(partitionings)
        parameterConvergencePlot(equation, integrators{i}, partitionings{j})
    end
end

%% Hyperviscosity Experiments
integrator = 'erk';
partitionings = {'hyperdiffusion-4', 'hyperdiffusion-6', 'hyperdiffusion-8'};
for j = 1 : length(partitionings)
    parameterConvergencePlot(equation, integrator, partitionings{j})
end