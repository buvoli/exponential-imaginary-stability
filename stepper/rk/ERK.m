function [ts, ys, tcpu, tccpu] = ERK(L, N, tspan, y0, Nt, options)
% ERK implements a generic exponential Runge Kutta integratorof the form:
%
%      Y_1 = N(t_n, y_n)
%      Y_i = d(1, i-1) * exp(h d(2, i-1) L) y_n + \sum_{j=1}^{i-1} \alpha_{ij}(hL) N(c_j, Y_j)    i = 2, ..., s
%
%      y_{n+1} = d(s) * exp(d(2, s) L) y_n + \sum_{j=1}^s \beta_{j}(hL) N(c_j, Y_j)
%
% where the functions
%
%      alpha_{ij}(L) = \sum_{k=1}^m A(1,k,j,i) * \varphi_k( A(2,k,j,i) * L )
%      beta_{j}(L)   = \sum_{k=1}^m b(1,k,j) * \varphi_k( b(2,k,j) * L)
%
% The information for the EXPRK A matrix is stored as a 4D array A(:,:,:,:) where A(:, :, j, i) specifies the linear
% combination of phi functions that that comprise the (i+1,j)th entry in the A matrix. In particular:
%
%          1. size(A, 3) = size(A, 4) is equal to (s-1) - first stage is explicit and not represented
%
%          2. size(A, 2) returns the highest-order phi-function used by the method
%
%          3. A(:, k, j, i) is a vector of dimension 2 that cooresponds to the term A(1, k, j, i) * phi_{k-1}(A(2, k, j, i) hL)        
%
% The b vector is stored a 3D array b(:,:,:) where b(:,:,i) contains the information for the ith entry of b.
%
%          1. size(b, 3) = s (the total number of stages)
%
%          2. b(:, k, j) is a vector of dimension 2 that cooresponds to the term b(1, k, j) * phi_{k-1}(b(2, k, j))
%
% The d vector is stored as a 2D array d(:) where d(:, i) cooresponds to the term d(1,j) exp( h * d(2, j) ) y_n
%
% The c vector is simple an vector of length s-1
%
% PARAMETERS
%   L       - vector, contains diagonal components of linear operator
%   N       - function, nonlinear operator
%   tspan   - integration bounds
%   y0      - initial condition
%   Nt      - number of timesteps
%   options - struct with fields: 
%               coeffGenerator     - function that returns exponential matrices A, b, c, and d               
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
[A, b, c, d] = options.coeffGenerator();
validateOptions(A, b, c, d);
tic;
[A_phi, b_phi, d_phi, A_ci] = initRKCoefficients(transpose(L(:)), h, A, b, d);
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
nstages = size(A, 3)  + 1;
Y = zeros(length(y0), nstages);
YN = zeros(length(y0), nstages);

tic;
for n = 1 : Nt
    % -- compute stage values ------------------------------------------------------------------------------------------
    Y(:,1) = y0;
    YN(:,1) = N(t0, y0, params);
    for i = 2 : nstages
        Y(:, i) = d_phi(:, i - 1) .* Y(:,1);
        for j = 1 : i - 1
            A_coeff_index = A_ci(j,i - 1);
            if(A_coeff_index ~= 0)
                Y(:, i) = Y(:, i) + A_phi(:, A_coeff_index) .* YN(:, j);
            end
        end
        YN(:, i) = N(t0 + h * c(i - 1), Y(:, i), params);
    end
    % -- compute output ------------------------------------------------------------------------------------------------
    y0 = d_phi(:, nstages) .* Y(:, 1);
    for i = 1 : nstages
        y0 = y0 + b_phi(:, i) .* YN(:, i);
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


function validateOptions(A, b, c, d)

s = size(A, 3) + 1; % number of total stages (including explicit stages)
k = size(A, 2); % number of phi functions used by method

if(( s == 2 && ~isequal(size(A), [2, k])) || ( s > 2 && ~isequal(size(A), [2, k, s - 1, s - 1])))
    error('invalid dimensions for A matrix')
end
if(( s == 1 && ~isequal(size(b), [2, k])) || (s > 1 && ~isequal(size(b), [2, k, s])))
    error('invalid dimensions for b matrix')
end
if(~isvector(c) || length(c) ~= (s - 1))
    error('invalid dimensions for c vector')
end
if(~isequal(size(d), [2 s]))
    error('invalid dimensions for c vector')
end

end


function [A_ci, nnz_A] = initAci(A)
% Form the Matrix A_coeffind. A_ci(i,j) is nonzero if exponential matrix coefficients in A(i,j) are nonempty (i.e. A(:,:,j,i) is non empy).
% The value cooresponds to the squential index of the nonzero coefficients. nnz_A is the total number of nonzero coefficients.

s   = size(A, 3) + 1; % num stages
A_ci  = zeros(s - 1, s - 1);
nnz_A = 0;

for i = 1 : (s - 1)
    for  j = 1 : (s - 1)
        if(any(A(1,:,j,i) ~= 0))
            nnz_A   = nnz_A + 1;
            A_ci(j,i) = nnz_A;
        end
    end
end

end


function [A_phi, b_phi, d_phi, A_ci] = initRKCoefficients(L, h, A, b, d)

nL      = length(L);
nphi    = size(A, 2);
nstages = size(A, 3) + 1;

[A_ci, A_nnz] = initAci(A);

% -- init A Tablaeu coefficients --------------------------------------------------------------------------------------- 
A_phi = zeros(nL, A_nnz);
for i = 1 : nstages - 1
    for j = 1 : i
        coeff_index = A_ci(j,i);
        if(coeff_index ~= 0) % A(i,j) is nonzero
            for k = 1 : nphi
                c1 = A(1, k, j, i);
                c2 = A(2, k, j, i);
                P = phi(h * c2 * L, nphi-1);
                A_phi(:, coeff_index) = A_phi(:, coeff_index) + h * c1 * P(k, :).';
                
            end
        end
    end
end

% -- init b vector coefficients ----------------------------------------------------------------------------------------
b_phi = zeros(nL, nstages);
for i = 1 : nstages
    for k = 1 : nphi
        c1 = b(1, k, i);
        c2 = b(2, k, i);
        if(c1 ~= 0) % b(i) contains \varphi_{k-1}(c2 L)
            P = phi(h * c2 * L, nphi - 1);
            b_phi(:, i) = b_phi(:, i) + h * c1 * P(k, :).';
        end
    end
end

% -- init d vector coefficients ----------------------------------------------------------------------------------------
d_phi = zeros(nL, nstages);
for i = 1 : nstages
    c1 = d(1, i);
    c2 = d(2, i);
    if(c1 ~= 0.0) % d(i) contains \exp(c2 L)
        P = phi(h * c2 * L, nphi-1);
        d_phi(:, i) = c1 * P(1, :).';
    end
end

end