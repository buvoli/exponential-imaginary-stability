addpath(genpath('../../stepper/common'));
addpath(genpath('../../stepper/rk'));

%% == load equation =======================================================
equation = 'zds';
wd = pwd();
cd(fullfile('../../stepper/equations/', equation))
init
cd(wd);

Nt = 2000;
h  = diff(tspan) / Nt;
ks = [0:pars.Nx/2 -pars.Nx/2+1:-1] * (2*pi/pars.Lx);

%% == Produce ERK Snapshot Plots ==========================================

[ts,ys] = ERK(LF(pars), NF, tspan, y0, Nt, struct('coeffGenerator', @ERK4, 'parameters', pars)); 
fh = solutionPlot(ts, xs, ys, filter, []);
exportFigure(fh, struct('SavePath', 'ERK-solution', 'Format', 'pdf'));

fh = spectrumPlot(ts, fftshift(ks), ys, sfilter, []);
exportFigure(fh, struct('SavePath', 'ERK-spectrum', 'Format', 'pdf'));

fprintf('rho(hL) = %f\n', max(h * abs(LF(pars))))

%% == Produce IMRK Snapshots ==============================================

[ts,ys] = IMRK(LF(pars), NF, tspan, y0, Nt, struct('coeffGenerator', @IMRK4, 'parameters', pars)); 
fh = solutionPlot(ts, xs, ys, filter, []);
exportFigure(fh, struct('SavePath', 'IMRK-solution', 'Format', 'pdf'));

fh = spectrumPlot(ts, fftshift(ks), ys, sfilter, []);
exportFigure(fh, struct('SavePath', 'IMRK-spectrum', 'Format', 'pdf'));

function fh = solutionPlot(ts, xs, ys, filter, title)

fh = figure();
surf(ts, xs, filter(ys)); shading interp;
axis tight;
zlim([0,2]);
view([-60 48]); grid off;  
zlabel('$\|u\|^2$', 'Interpreter', 'latex');
ylabel('$x$', 'Interpreter', 'latex'); 
xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18, 'FontName', 'Minion Pro');

if(~isempty(title))
    title(title);
end

end

function fh = spectrumPlot(ts, ks, ys, filter, title)

fh = figure();
surf(ts, ks, filter(ys)); shading interp;
axis tight;
view([-60 48]); grid off;  
zlabel('$\ln(|a_k|)$', 'Interpreter', 'latex');
ylabel('$k$', 'Interpreter', 'latex'); xlabel('t', 'Interpreter', 'latex');
set(gca, 'FontSize', 18, 'FontName', 'Minion Pro');


if(~isempty(title))
    title(title);
end

end