function [R] = rERK1(z1, z2_vec)
%ampERK1 stability function for ERK1
% y_{n+1} = y_{n} + h \varphi_1(h L) F(y_{n})
%
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
P = phi(z1, 1);
A = 0;
d = [1; P(1)];
b = P(2);
end