addpath(genpath('../../stepper/common'));
addpath(genpath('../../stepper/rk'));

%% == load equation =======================================================
equation = 'kdv';
wd = pwd();
cd(fullfile('../../stepper/equations/', equation))
init
cd(wd);

%% == Produce ERK Snapshot Plots ==========================================
Nt = ceil(2000 * tspan(end)); pars.rho = pi/128;
[ts, ys, tcpu] = ERK(LF(pars), NF, tspan, y0, Nt, struct('coeffGenerator', @ERK4, 'parameters', pars, 'max_ts_to_store', 2000));

%% == Produce ERK Snapshot Plots ==========================================
fh = figure();
surf(xs, ts, transpose(filter(ys))); 
shading interp; view([90 90]); 
axis([0 pars.Lx, tspan]); xticks([0 pars.Lx / 2 pars.Lx]);
ylabel('Time $t$', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
colorbar('northoutside')
set(gca, 'FontName', 'Minion Pro', 'FontSize', 14);

exportFigure(fh, struct('SavePath', 'kdv-solution', 'Format', 'pdf', 'PaperPosition',   [0 0 24 7]))