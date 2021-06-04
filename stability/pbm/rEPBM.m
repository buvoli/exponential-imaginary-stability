function [amp] = ampETDPBM(z1, z2_vec, options)
%ampETDPBM amplification factor amp = \rho(M(z_1, z_2)) for Parallel Exponential Polynomial Block Method
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
%   method_handle - (function_handle) handle to desired method
%   options - struct with fields:
%               z                   - nx1 array which contains pbm nodes
%               b                   - nx1 array which contains left endpoints
%               alpha               - extrapolation parameter
%               S                   - set S for PMFOmS (z(S) will not be used to form L_N
%               m                   - number of correction iterations per step
% RETURNS
%   amp    - amplification factor

% -- set default params ------------------------------------------------------------------------------------------------
default_value_cell = {
    {'z',     -cos(pi*((0:3)')/3)}
    {'b',     -1 * ones(3, 1)}
    {'alpha', 1}
    {'kappa', 0}
    {'mS',    []}
};
options = setDefaultOptions(options, default_value_cell);

% -- varify paramters are valid ----------------------------------------------------------------------------------------
if(length(options.z) ~= length(options.b))
    error('length(b) is not equal to length(a)');
end
if(options.z(1) ~= -1)
    error('z(1) must be negative one');
end

% -- Numerical Parameters ----------------------------------------------------------------------------------------------
z = options.z;
q = length(z);
b = options.b;
kappa = options.kappa;
LN_active_inds = setdiff(1 : q, options.mS);
alpha = options.alpha;

% === Initialize Coefficients & Block Matrices =========================================================================
A_method = zeros(q);
A_itr    = zeros(q); A_itr(:, 1) = 1;
for i = 1 : q
    A_method(i, :) = weights(b(i) / alpha, z / alpha, 0);
end    
[W_method, P0_method] = initW(z1, z(LN_active_inds) / alpha, b / alpha, z/alpha+1);
[W_itr,    P0_iter]    = initW(z1, z(LN_active_inds) / alpha, z(1) / alpha * ones(q, 1), z / alpha);    
    
% On Dalquist problem method reduces to y_{n+1} = (A(z1) + z2 B(z1)) * y_n
A_z1 = P0_method(:) .* A_method; % equivalent to diag(z1 * ones(q,1)) * A
B_z1 = zeros(q);
B_z1(:, LN_active_inds) = W_method;

iA_z1 = P0_iter(:) .* A_itr; 
iB_z1 = zeros(q);
iB_z1(:, LN_active_inds) = (W_itr);
  
num_z2 = length(z2_vec);
amp    = zeros(num_z2, 1);
for i = 1 : num_z2    
    z2 = z2_vec(i);
    if(kappa > 0)
        M_z1 = (iA_z1 + z2 * iB_z1)^kappa * (A_z1 + z2 * B_z1);
    else
        M_z1 = (A_z1 + z2 * B_z1);
    end
    amp(i) = max(abs(eig(M_z1)));
end

end

function [W,P0]=initW(z1,tau,a,b)
%INITW Initializes ETDMP coefficients for scaler L using weights function by Fornberg.
% Input Parameters
%   z1   - (scalar) lambda_1
%   tau - (array) quadrature points used for MP method
%   a   - (array) left integration bounds
%   b   - (array) right integration bounds
% Output Parameters
%   W   - (array) array of dimensions length(L) x n x n.
%         W(:, :,j) contains the coefficients for computing the jth output
eta = b - a;
q   = length(tau);
n   = length(a);
W   = zeros(n, q);
P0  = zeros(q, 1);             % stores phi_0
for i=1:n
    p = phi(eta(i)*z1,q);                % phi functions 0 to n
    stau = (tau - a(i))/eta(i);              % scaled quadrature ppints
    if(eta(i) ~= 0)
        A = weights(0,stau,q-1);              % finite difference weights
    else
        A = zeros(q,q);
    end
    W(i, :) = eta(i)*A*p(2:end,:);     % store ith row of W matrix
    P0(i)    = p(1,:);                               % store phi_0
end
end