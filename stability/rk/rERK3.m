function [R] = rERK3(z1, z2_vec)
%rERK3 stability function for Eqn (23)-(25) from S. M. Cox and P. C. Matthews, "Exponential time differencing for stiff systems." 2002 
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
P = phi(z1, 4);
P_12 = phi(z1/2, 2);

d      = zeros(4,1);
d(1)   = 1;                             % 1
d(2)   = P_12(1);                       % exp(h/2 L)
d(3:4) = P(1);                          % exp(h L)

A = zeros(3);
A(1,1) = 1;
A(2,1) = (1 / 2) * P_12(2);             % (1/2) \varphi_1(hL/2)
A(3,1) = (-1) * P(2);                   % -1 \varphi_1(hL) 
A(3,2) = (2) * P(2);                    % 2  \varphi_1(hL)

b    = zeros(3,1);
b(1) = (P(2) - 3 * P(3) + 4 * P(4));    % (\varphi_1(hL) - 3 \varphi_2(hL) + 4 \varphi_3(hL)) * F(Y1)
b(2) = (4 * P(3) - 8 * P(4));           % (4 \varphi_2(hL) - 8 \varphi_3(hL)) * F(Y2)
b(3) = (-1 * P(3) + 4 * P(4));          % (- \varphi_2(hL) + 4 \varphi_3(hL)) * F(Y3) 
end