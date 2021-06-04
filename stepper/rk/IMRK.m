function [ts, ys, tcpu, tccpu] = IMRK(L, N, tspan, y0, Nt, options)
% IMRK implements a generic IMEX Runge-Kutta integrator of the form:
%
%      Y_i = y_n + \sum_{j=1}^{i-1} A(i,j) F(c_j, Y_j)    i = 1, ..., s
%      y_{n+1} = d(i) y_n + \sum_{j=1}^s b(i) F(c_j, Y_j)
%
% The IMRK matrices are produced by calling
% [A_im, b_im, A_ex, b_ex, c] = options.coeffGenerator();
%
% PARAMETERS
%   L       - vector, contains diagonal components of linear operator
%   N       - function, nonlinear operator
%   tspan   - integration bounds
%   y0      - initial condition
%   Nt      - number of timesteps
%   options - struct with fields: 
%               coeffGenerator     - function that returns implicit and explicit matrices A, b, and c               
%               max_ts_to_store    - max solution values to store
%               problem_parameters - struct which is passed to L and N Functions

% === START Load Options ===============================================================================================
default_value_cell = {
    {'max_ts_to_store', max(2,min(5000000/length(y0),1000))}
};
options = setDefaultOptions(options, default_value_cell);

% -- varify paramters are valid ----------------------------------------------------------------------------------------
if(~isfield(options, 'parameters'))
    error('options.parameters was not provided');
end

max_ts_to_store = max(2, options.max_ts_to_store);
params  = options.parameters;
% === END Load Options =================================================================================================

% Numerical Parameters
h = (tspan(end)-tspan(1))/Nt;

% === START Initialize Coefficients ====================================================================================
[A_im, b_im, A_ex, b_ex, c] = options.coeffGenerator();
validateOptions(A_im, b_im, A_ex, b_ex, c);
tic;
tccpu = toc;
% === END Initialize Coefficients ======================================================================================

% Data Parameters
skip_rate = ceil(Nt/max_ts_to_store);
ys   = zeros(length(y0),ceil(Nt/skip_rate)+1);
ts   = zeros(size(ys,2),1);
ys(:,1) = y0; save_count = 2;
ts(1)   = tspan(1); t0 = tspan(1);
y0      = reshape(y0, length(y0), 1);

% ---- stepping code ---------------------------------------------------------------------------------------------------
nstages = size(A_im, 1);
Y = zeros(length(y0), nstages);
YL = zeros(length(y0), nstages);
YN = zeros(length(y0), nstages);
L = L(:); % turn into row vector

tic;
for n = 1 : Nt
    % -- compute stage values ------------------------------------------------------------------------------------------
    for i = 1 : nstages
        Y(:, i) = y0;
        for j = 1 : (i - 1)
            Y(:, i) = Y(:, i) + (h * A_im(i,j)) * (YL(:, j)) + (h * A_ex(i,j)) * YN(:, j);
        end
        Y(:, i) = Y(:, i) ./ ( 1 - (h * A_im(i,i)) * L);
        % -- evaluate stage derivatives
        YL(:,i) = L .* Y(:, i);
        YN(:,i) = N(t0 + h * c(i), Y(:,i), params);
    end
    % -- compute output ------------------------------------------------------------------------------------------------
    for i = 1 : nstages
        y0 = y0 + (h * b_im(i)) .* YL(:, i) + (h * b_ex(i)) .* YN(:, i);
    end
    t0 = t0 + h;
    % Save Data
    if(mod(n,skip_rate) == 0 || n==Nt)
        ys(:,save_count) = y0;
        ts(save_count,:) = t0;
        save_count = save_count + 1;
    end    
end
tcpu = toc;

end


function validateOptions(A_im, b_im, A_ex, b_ex, c)

s = size(A_im, 1); % number of total stages (including explicit stages)
dims = [size(A_im, 2) size(A_ex, 1) size(A_ex, 2) length(b_im) length(b_ex) length(c)];

if(~all(dims == s))
    error('invalid coefficient matrices: A_im, A_ex should be sxs while b_im, b_ex, and c should have length s')
end

end
