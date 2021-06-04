function [R] = rERK2(z1, z2_vec)
%rERK2 stability function for Eqn (22) from S. M. Cox and P. C. Matthews, "Exponential time differencing for stiff systems." 2002
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
% RETURNS
%   R    - value of stability function

[d, A, b] = initETDCoefficients(z1);
R = rERK(z2_vec, d, A, b);
end

function [d, A, b] = initETDCoefficients(z1)
%INITETDRK4 Initializes coefficients for ERK1 for scaler Lambda
P = phi(z1, 3);
d = [1; P(1); P(1)];
A = [ 0 0; P(2) 0 ];
b = [P(2) - P(3); P(3)];
end