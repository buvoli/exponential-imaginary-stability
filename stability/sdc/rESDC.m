function amp = rESDC(z1, z2_vec, params)
%ampESDC amplification factor amp = \rho(M(z_1, z_2)) for exponential sdc
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2_vec (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
%   params - struct with fields:
%               z                   - nx1 array which contains sdc nodes
%               m                   - nx1 array which contains left endpoints
% RETURNS
%   amp - amplification factor
    
[I, P01] = initW(z1, 1, params.tau);
amp = ampF(length(params.tau), params.m, z2_vec, P01, I);

end

function psi = ampF(n, m, z2, P01, I)
%AMP_ETD Computes amplification factor (labeled \psi(r,z) in paper) for ETD 
% spectral deferred correction scheme assuming n chebyshyshev quadrature 
% points and m corrections sweeps.
%PARAMETERS
% n  - number of quadrature points
% m  - number of correction sweeps
% z2 - explicit component
% I  - ETD integration matrix

z2 = z2(:).';

% Euler
psi = ones(n,length(z2)); 
for i=1:n-1
    psi(i+1,:) = (P01(i,1) + z2*P01(i,2)).*psi(i,:);   
end

% Correction Sweeps
for k=1:m
    psi_old = psi;
    for i=1:n-1
        psi(i+1,:) = P01(i,1)*psi(i,:) + P01(i,2)*z2.*(psi(i,:) - psi_old(i,:)) + z2.*(I(i,:)*psi_old);
    end
end
psi = psi(end, :);
end

function [W, P01] = initW(L, h, tau)
%INITW_SCALER Initializes ETDSDC W functions for scaler/vector L using weights function by Fornberg.
% Input Parameters
%   L   - (array)  L = [\Lambda_1, \Lambda_2, ..., \Lambda_N]
%   h   - (double) timestep
%   tau - (array)  normalized quadrature points
% Output Parameters
%   W   - (array) array of dimensions (n-1) x n x length(L). 
%         W(:,:,j) contains integration matrix cooresponding to L(j)
%   P01 - (array) array containg phi_0 and phi_1 needed for ETD Euler method.

eta = tau(2:end) - tau(1:end-1);
n   = length(tau);
W   = zeros(n-1,n,length(L));             % stores integration matrices W
P01 = zeros(n-1,2,length(L));             % stores phi_0 and phi_1 for ETD Euler
for i=1:n-1
    q = (tau - tau(i))/eta(i);            % scaled quadrature ppints
    a = weights(0,q,n-1);                 % finite difference matrix
    p = phi(eta(i)*h*L,n);         % phi functions 0 to n          
    W(i,:,:) = h*eta(i)*a*p(2:end,:);     % store ith row of W matrix
    P01(i,1,:) = p(1,:);                  % store phi_0 and phi_1
    P01(i,2,:) = h*eta(i)*p(2,:);         % store phi_0 and phi_1
end
end