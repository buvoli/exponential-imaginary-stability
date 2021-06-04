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

%%

ax = [-120 0 -1500 1500];

f = figure(1);
rhos = pi ./ ( 2 .^ [Inf 11 9 7 5] );
for i = 1 : length(rhos)
    L = LF(mergeStructs(pars, struct('rho', rhos(i)))); [~, inds]=sort(imag(L)); L = L(inds);
    plot(real(L), imag(L), '.-', 'MarkerSize', 14, 'LineWidth', 1.5); hold on;
end
hold off;
axis(ax);
xlabel('Re($z$)', 'interpreter', 'Latex'); ylabel('Im($z$)', 'interpreter', 'Latex');
legend(arrayfun(@rhoLabel, rhos, 'UniformOutput', false),'Location', 'northoutside', 'Orientation', 'Horizontal', 'Interpreter', 'Latex'); legend box off;
set(gca, 'FontSize', 12, 'FontName', 'Minion Pro');
exportFigure(f, struct('SavePath', ['output/spectrum/spectrum-', equation, '-rotate'], 'Format', 'pdf', 'PaperPosition', [0 0 20 6], 'Renderer', 'opengl'));

f = figure(2);
rhos = pi ./ ( 2 .^ [Inf 8 6 4 2] );
for i = 1 : length(rhos)
    L = LF(mergeStructs(pars, struct('rhouxx', rhos(i)))); [~, inds]=sort(imag(L)); L = L(inds);
    plot(real(L), imag(L), '.-', 'MarkerSize', 14, 'LineWidth', 1.5); hold on;
end
hold off;
axis(ax);
xlabel('Re($z$)', 'interpreter', 'Latex'); ylabel('Im($z$)', 'interpreter', 'Latex');
legend(arrayfun(@rhoLabel, rhos, 'UniformOutput', false),'Location', 'northoutside', 'Orientation', 'Horizontal', 'Interpreter', 'Latex'); legend box off;
set(gca, 'FontSize', 12, 'FontName', 'Minion Pro');
exportFigure(f, struct('SavePath', ['output/spectrum/spectrum-', equation, '-underrotate'], 'Format', 'pdf', 'PaperPosition', [0 0 20 6], 'Renderer', 'opengl'));

f = figure(3);
epsilons = [0 1 2 4 8];
for i = 1 : length(epsilons)
    L = LF(mergeStructs(pars, struct('epsilon', epsilons(i)))); [~, inds]=sort(imag(L)); L = L(inds); L = L(L ~= 0);
    plot(real(L), imag(L), '.-', 'MarkerSize', 14, 'LineWidth', 1.5); hold on;
end
hold off;
axis(ax);
xlabel('Re(z)'); ylabel('Im(z)');
legend(arrayfun(@epsilonLabel, epsilons, 'UniformOutput', false),'Location', 'northoutside', 'Orientation', 'Horizontal', 'Interpreter', 'Latex'); legend box off;
set(gca, 'FontSize', 12, 'FontName', 'Minion Pro');
exportFigure(f, struct('SavePath', ['output/spectrum/spectrum-', equation, '-translate'], 'Format', 'pdf', 'PaperPosition', [0 0 20 6], 'Renderer', 'opengl'));

function l = rhoLabel(rho)
    if(rho == 0)
        l = ['$\rho = 0 \hspace{1em}$'];
    else
        [~, b] = rat(rho / pi);
        l = ['$\rho = \pi / ', num2str(b), '\hspace{1em}$'];
    end
end

function l = epsilonLabel(epsilon)
    l = ['$ \epsilon = ', num2str(epsilon), '\hspace{1em}$'];
end