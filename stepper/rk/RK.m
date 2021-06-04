function [t0, y0, ts, ys] = RK(f, tspan, y0, Nt, options)
% RK implements a generic explicit Runge Kutta integrator of the form:
%
%      Y_1 = F(t_n, y_n)
%      Y_i = y_n + \sum_{j=1}^{i-1} h a_{ij} F(c_j, Y_j)    i = 2, ..., s
%
%      y_{n+1} = y_n + \sum_{j=1}^s h b_j F(c_j, Y_j)
%
% The RK matrices are produced by calling
%   [A, b, c] = options.coeffGenerator();
%
% PARAMETERS
%   f       - rhs function
%   tspan   - integration bounds
%   y0      - initial condition
%   Nt      - number of timesteps
%   options - struct with fields: 
%               coeffGenerator     - function that returns RK matrices A, b, and c               
%               max_ts_to_store    - max solution values to store
%               problem_parameters - struct which is passed to L and N Functions
 
t0 = tspan(1);
h  = diff(tspan) / Nt;
pars = options.parameters;

if(nargout >= 3)
    store_data = true;
    if(isfield(options, 'max_ts_to_store'))
        max_ts_to_store = options.max_ts_to_store;
    else
        max_ts_to_store = max(2,min(5000000/length(y0),1000));
    end
    skip_rate = ceil(Nt/max_ts_to_store);
    ys = zeros(length(y0),ceil(Nt/skip_rate)+1);
    ts = zeros(size(ys,2),1);
    ys(:,1) = y0; save_count = 2;
    ts(1) = tspan(1); t0 = tspan(1);
else
    store_data = false;
end

% === START Initialize Coefficients ====================================================================================
[A, b, c] = options.coeffGenerator();
validateOptions(A, b, c);
% === END Initialize Coefficients ======================================================================================

nstages = size(A, 1) + 1;
Y = zeros(length(y0), nstages);
YF = zeros(length(y0), nstages);

tic;
for n = 1 : Nt
    % -- compute stage values ------------------------------------------------------------------------------------------
    Y(:,1) = y0;
    YF(:,1) = f(t0, y0, pars);
    for i = 2 : nstages
       Y(:, i) = Y(:, 1);
        for j = 1 : i - 1
            Y(:, i) = Y(:, i) + h * A(i - 1, j) .* YF(:, j);
        end        
        YF(:, i) = f(t0 + h * c(i - 1), Y(:, i), pars);
    end
    % -- compute output ------------------------------------------------------------------------------------------------
    for i = 1 : nstages
        y0 = y0 + h * b(i) .* YF(:, i);
    end
    t0 = t0 + h;
    
    % -- Save Data -----------------------------------------------------------------------------------------------------
    if(store_data && (mod(i,skip_rate) == 0 || i==Nt))
        ys(:,save_count) = y0;
        ts(save_count,:) = t0;
        save_count = save_count + 1;
    end 
        
end
 
end

function validateOptions(A, b, c)

s = size(A, 1) + 1; % number of total stages (including explicit stages)
dims = [size(A, 2) length(b) length(c)];

if(~isequal(dims, [s-1, s, s-1]))
    error('invalid coefficient matrices: matrix A should be (s-1)x(s-1), and vectors b,c should have length s, and s-1.')
end

end