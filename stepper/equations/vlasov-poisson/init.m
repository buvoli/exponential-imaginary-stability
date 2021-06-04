% Spatial Descretization
Lx = 20*pi; Nx = 135; % 256 * 2; %135; %135; % 256 * 2; %135;
Lv = 16;    Nv = 256; % 256 * 2; %256; % 256 * 2;
xs = linspace(0,Lx,Nx+1); xs(end) = [];
vs = linspace(-Lv/2,Lv/2,Nv+1); vs(end) = [];
[X, V] = meshgrid(xs,vs);

% Temporal Discretization
tspan = [0 50];

% Parameters & Differentiation Matrices
pars = struct('Lx', Lx, 'Nx', Nx, 'Lv', Lv, 'Nv', Nv, 'vdes', 'W5', 'epsilon', 0, 'rho', 0/64, 'rhouxx', 0/32, 'antialias', false);
[DX, DV, IDX, DXX] = DM(pars); pars.DX = DX; pars.DV = DV; pars.IDX = IDX; pars.DXX = DXX; pars.E = @ETrap;

% initial conditions
y0 = ( 0.9 / sqrt(2 * pi) * exp(-V.^2  / 2) + 0.2 / sqrt(2 * pi) * exp(-2 * (V - 4.5).^2)) .* ( 1 + 0.04 * cos(0.3 * X));
y0 = fft(y0, [], 2); % fft in x direction only
y0 = reshape(y0,Nx*Nv,1);

% set L and N
LF = @L;
NF = @N;

% Plot/Error Variables
to_physical   = @(x) real(ifft(reshape(x,pars.Nv,pars.Nx), [], 2));
to_physical_v = @(x) reshape(to_physical(x),[pars.Nv * pars.Nx, 1]);
to_fourier    = @(x) fft(reshape(x,pars.Nv,pars.Nx), [], 1);

% Performance Experiment Settings
Nt_reference = 60000;
Nts = round(logspace(2,4,12));

error_norm = @(x,y) max(abs(to_physical_v(x-y)))./max(abs(to_physical_v(y)));

%% Example Run
% mesh(xs, vs, to_physical(y0))
%Nt = 800; [ts, ys] = IMRK(LF(pars), NF, tspan, y0, Nt, struct('coeffGenerator', @IMRK4, 'parameters', pars));
%Nt = 800; [ts, ys, tcpu] = ERK(LF(pars), NF, tspan, y0, Nt, struct('coeffGenerator', @ERK4, 'parameters', pars));
%surf(xs, vs, to_physical(ys(:,end))); view([90 90]); axis tight; shading interp; colorbar; caxis([0 0.4]);


%% Energy Drift Test
% ns = size(ys,2);
% en = zeros(ns, 2);
% 
% for i = 1 : ns
%     [en(i,1), en(i,2)] = HTrap(ys(:,i), pars);
% end
% en(:,1) = (en(:,1) - en(1,1)) / en(1,1);
% 
% figure(); plot(ts, en(:,1)); title('Relative energy Drift');
% figure(); plot(ts, en(:, 2)); title('energy \|E\|^2');