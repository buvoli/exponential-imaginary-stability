function [R] = rERK4(z1, z2_vec)
%rERK4 stability function for Eqn (51) from S. Krogstad. "Generalized integrating factor methods for stiff PDEs", 2005
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

P = phi(z1, 5);
P_12 = phi(z1 / 2, 2);

d = zeros(5, 1);
d(1)   = 1;                             % 1
d(2:3) = P_12(1);                       % exp(h/2 L)
d(4:5) = P(1);                          % exp(h L)

A = zeros(4);
A(2,1) = (1/2) * P_12(2);               % (1/2) \varphi_1(h/2 L)
A(3,1) = (1/2) * P_12(2) - P_12(3);     % (1/2) \varphi_1(h/2 L) - \varphi_2(h/2 L)
A(3,2) = P_12(3);                       % \varphi_2(h/2 L)
A(4,1) = P(2) - 2 * P(3);               % \varphi_1(h L) - 2 \varphi_2(h L)
A(4,3) = 2 * P(3);                      % 2 \varphi_2(h L)

b = zeros(4, 1);
b(1) = P(2) - 3 * P(3) + 4 * P(4);      % (\varphi_1(hL) - 3 \varphi_2(hL) + 4 \varphi_3(hL) * F(y_n)
b(2) = 2 * P(3) - 4 * P(4);             % (2 \varphi_2(hL) - 4 \varphi_3(hL)) * F(Y_1)
b(3) = 2 * P(3) - 4 * P(4);             % (2 \varphi_2(hL) - 4 \varphi_3(hL)) * F(Y_2)
b(4) = -1 * P(3) + 4 * P(4);            % (-1 \varphi_2(hL) - 4 \varphi_3(hL)) * F(Y_3)

end