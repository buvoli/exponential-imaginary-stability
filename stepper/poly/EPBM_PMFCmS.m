function [ts, ys, tcpu, tcpuc] = EPBM_PMFCmS(L, N, tspan, yi, Nt, options)
%EPBM Parallel Exponential Polynomial Block Method (PMFO and PMFOmS)
% PARAMETERS
%   L       - vector, contains diagonal components of linear operator
%   N       - function, nonlinear operator
%   tspan   - integration bounds
%   yi      - initial condition
%   options - struct with fields:
%               max_ts_to_store     - max solution values to store
%               problem_parameters  - struct which is passed to L and N Functions
%               z                   - nx1 array which contains pbm nodes
%               b                   - nx1 array which contains left endpoints
%               alpha               - extrapolation parameter
%               S                   - set S for PMFOmS (z(S) will not be used to form L_N)
% RETURNS
%   ys      - solution values
%   ts      - times cooresponging to solution values
%   tcpu    - seconds to compute solution
%   tcpuc   - second to initialize coefficients

% === START Load Options ===============================================================================================
default_value_cell = {
    {'z',     -cos(pi*((0:3)')/3)}
    {'b',     -1 * ones(3, 1)}
    {'alpha', 1}
    {'kappa',     0}
    {'mS',    []}
    {'max_ts_to_store', max(2,min(5000000/length(yi),1000))}
};
options = setDefaultOptions(options, default_value_cell);

% -- varify paramters are valid ----------------------------------------------------------------------------------------
if(length(options.z) ~= length(options.b))
    error('length(b) is not equal to length(a)');
end
if(options.z(1) ~= -1)
    error('z(1) must be negative one');
end
if(~isfield(options, 'parameters'))
    error('options.parameters was not provided');
end

max_ts_to_store = max(2, options.max_ts_to_store);
params  = options.parameters;
% === END Load Options =================================================================================================

% -- Numerical Parameters ----------------------------------------------------------------------------------------------
h  = (tspan(end)-tspan(1))/Nt;
z = options.z;
q = length(z);
b = options.b;
kappa = options.kappa;
LN_active_inds = setdiff(1 : q, options.mS);
alpha = options.alpha;

Nx = length(yi);
y_n   = repmat(yi(:), 1, q); % constant initial guess
y_np1 = zeros(size(y_n));
t_n   = tspan(1);

% Data Parameters
if(nargin == 7)
    max_ts_to_store = max(2,min(5000000/length(yi),1000));
end
skip_rate = ceil(Nt/ max_ts_to_store);
ys      = zeros(Nx,ceil(Nt/skip_rate)+1);
ts      = zeros(size(ys,2),1);
ys(:,1) = yi; save_count = 2;
ts(1)   = tspan(1);

% === START Initialize Coefficients ====================================================================================
tic;
AT = zeros(q);
for i = 1 : q
    AT(:, i) = weights(b(i) / alpha, z / alpha, 0);
end
[W_mthd,P0_mthd] = initW(transpose(L(:)), h, z(LN_active_inds) / alpha, b / alpha, z/alpha+1);
[W_itr, P0_itr]  = initW(transpose(L(:)), h, z(LN_active_inds) / alpha, z(1) / alpha * ones(q, 1), z / alpha);
tcpuc = toc;
% === END Initialize Coefficients ======================================================================================

% === START Initial Condition ==========================================================================================
for k = 1 : q + kappa   % Run Polynomial Iterator
    N_n = N(t_n + h*((z+1)/alpha), y_n(:, LN_active_inds), params);
    for j = 1 : q
        y_np1(:, j) = P0_itr(:, j) .* y_n(:, 1) + sum(W_itr(:,:,j) .* N_n, 2);
    end
    y_n = y_np1;
end
% ==== End Initial Condition ===========================================================================================

tic; % TIME STEPPING LOOP
for i = 1 : Nt
    N_n = N(t_n + h*((z+1)/alpha),y_n(:, LN_active_inds),params);
    for j = 1 : q
        y_np1(:, j) = P0_mthd(:, j) .* (y_n * AT(:, j)) + sum(W_mthd(:,:,j) .* N_n, 2);
    end
    y_n = y_np1;
    for k = 1 : kappa
        N_n = N(t_n + h*((z+1)/alpha) + h, y_n(:, LN_active_inds), params);
        for j = 1 : q
            y_np1(:, j) = P0_itr(:, j) .* y_n(:, 1) + sum(W_itr(:,:,j) .* N_n, 2);
        end
        y_n = y_np1;
    end
    
    t_n = t_n + h;
    % Save Data
    if(mod(i,skip_rate) == 0 || i == Nt)
        ys(:,save_count) = y_n(:, 1);
        ts(save_count,:) = t_n;
        save_count = save_count + 1;
    end
end
tcpu=toc;
end


function [W,P0]=initW(L,h,tau,a,b)
%INITW Initializes ETDMP coefficients for scaler L using weights function by Fornberg.
% Input Parameters
%   L   - (array) of \Lambda values L = [\Lambda_1, \Lambda_2, ..., \Lambda_N]
%   tau - (array) quadrature points used for MP method
%   a   - (array) left integration bounds
%   b   - (array) right integration bounds
% Output Parameters
%   W   - (array) array of dimensions length(L) x n x n.
%         W(:, :,j) contains the coefficients for computing the jth output
eta = b - a;
q   = length(tau);
n   = length(a);
W   = zeros(length(L), q, n);
P0  = zeros(length(L), q);             % stores phi_0
for i=1:n
    p = phi(h*eta(i)*L,q);                % phi functions 0 to n
    stau = (tau - a(i))/eta(i);              % scaled quadrature ppints
    if(eta(i) ~= 0)
        A = weights(0,stau,q-1);              % finite difference weights
    else
        A = zeros(q,q);
    end
    W(:, : ,i) = transpose(h*eta(i)*A*p(2:end,:));     % store ith row of W matrix
    P0(:,i)    = p(1,:);                               % store phi_0
end
end