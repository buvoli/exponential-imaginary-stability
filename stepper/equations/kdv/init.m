% Spatial Descretization
Lx = 2; Nx = 2^8;
xs = linspace(0,Lx,Nx+1)'; xs(end) = [];

% Temporal Discretization
sf = 1.6;
tspan = [0 100] * sf;

% Paremeters
pars = struct('Lx',Lx,'Nx',Nx,'delta',.022, 'rho', 0, 'rhouxx', 0, 'epsilon', 0, 'epsilon_tr', 0/10, 'epsilon_rt', 0/10000, 'antialias', true); % epsilon_tr dmp stabilizes ERK (approx 30,000 to 120,000 timesteps)

% Initial conditions
y0 = fft(cos(pi*xs));
if(pars.antialias)
    y0 = antialias(y0);
end

% set L and N
LF = @L;
NF = @N;

% Experiment Settings
Nt_reference = 750000 * sf;
Nt_run = 35000 * sf;
Nts = round(logspace(2,6,12) * sf);

% Solution filter
filter = @(x) abs(ifft(x));
sfilter = @(x) fftshift(log(abs(x)),1);
error_norm = @(x,y) max(abs(real(ifft(x-y)))./abs(real(ifft(y))));