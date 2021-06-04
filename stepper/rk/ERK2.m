function [A, b, c, d] = ERK2()
%ERK2 Eqn (22) from S. M. Cox and P. C. Matthews, "Exponential time differencing for stiff systems." 2002 

s = 2; % number of stages (including explicit first stage)
m = 3; % number of phi functions (including \varphi_0)

A = zeros(2, m, s - 1, s - 1);
b = zeros(2, m, s);
d = zeros(2, s);
c = zeros(s - 1, 1); 

% -- A Matrix ------------------------------------------------------------------------
A(:, 2, 1, 1) = [1.0,  1.0]; % \varphi_1(hL)

% -- b Vector ------------------------------------------------------------------------
% (\varphi_1(hL) - \varphi_2(hL)) * F(Y1)
b(:, 3, 1)    = [-1.0, 1.0];
b(:, 2, 1)    = [1.0,  1.0]; 

% \varphi_2(h L) * F(Y2)
b(:, 3, 2)    = [1.0,  1.0];

% -- d Vector ------------------------------------------------------------------------
d(:, 1)       = [1.0,  1.0]; % exp(h L) y_n
d(:, 2)       = [1.0,  1.0]; % exp(h L) y_n

% -- c Vector ------------------------------------------------------------------------
c(:)          = [1.0];


end