% Spatial Descretization
Lx = 8 * pi; Nx = 2 ^ 7;
xs = linspace(-Lx/2, Lx/2, Nx+1); xs(end) = [];

% Temporal Discretization
tspan = [0 40];

% initial conditions
mn = 3;
y0 = fft(1 + (1/100) * exp(2*pi*1i*mn*xs/Lx));
pars = struct('Lx', Lx, 'Nx', Nx, 'antialias', true, 'rho', 0, 'rhouxx', 0, 'epsilon', 0, 'method_stepsize', 0, 'method_order', 5, 'hyperv_order', 0, 'hyperv_coeff', 0); % epsilon 10 is required for ERK, or delta = 1/100
if(pars.antialias)
    y0 = antialias(y0);
end

% set L and N
LF = @L;
NF = @N;

% Performance Experiment Settings
Nt_reference = 200000;
Nts = round(logspace(2,5,12));

% Solution filter
filter = @(x) abs(ifft(x));
sfilter = @(x) fftshift(log(abs(x)),1);
error_norm = @(x,y) max(abs(real(ifft(x-y)))./abs(real(ifft(y))));