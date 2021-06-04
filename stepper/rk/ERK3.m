function [A, b, c, d] = ERK3()
%ERK3 Eqn (23)-(25) from S. M. Cox and P. C. Matthews, "Exponential time
% differencing for stiff systems." 2002

s = 3; % number of stages (including explicit first stage)
m = 4; % number of phi functions (including \varphi_0)

A = zeros(2, m, s - 1, s - 1);
b = zeros(2, m, s);
d = zeros(2, s);
c = zeros(s - 1, 1);

% -- A Matrix ------------------------------------------------------------------------
A(:, 2, 1, 1) = [0.5,  0.5]; % (1/2) \varphi_1(hL/2)
A(:, 2, 1, 2) = [-1.0, 1.0]; % -1 \varphi_1(hL)
A(:, 2, 2, 2) = [2.0,  1.0]; % 2  \varphi_1(hL)

% -- b Vector ------------------------------------------------------------------------
% (\varphi_1(hL) - 3 \varphi_2(hL) + 4 \varphi_3(hL)) * F(Y1)
b(:, 2, 1)    = [1.0, 1.0];
b(:, 3, 1)    = [-3.0, 1.0];
b(:, 4, 1)    = [4.0, 1.0];
% (4 \varphi_2(hL) - 8 \varphi_3(hL)) * F(Y2)
b(:, 3, 2)    = [4.0, 1.0];
b(:, 4, 2)    = [-8.0, 1.0];
% (- \varphi_2(hL) + 4 \varphi_3(hL)) * F(Y3)
b(:, 3, 3)    = [-1.0, 1.0];
b(:, 4, 3)    = [4.0, 1.0];

% -- d Vector ------------------------------------------------------------------------
d(:, 1)       = [1.0,  0.5]; % exp(h L / 2) y_n
d(:, 2)       = [1.0,  1.0]; % exp(h L) y_n
d(:, 3)       = [1.0,  1.0]; % exp(h L) y_n

% -- c Vector ------------------------------------------------------------------------
c(:)          = [0.5, 1.0];

end
