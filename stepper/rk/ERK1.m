function [A, b, c, d] = ERK1()
%ERK1 Exponential Euler

s = 1; % number of stages (including explicit first stage)
m = 2; % number of phi functions (including \varphi_0)

A = zeros(2, m, s - 1, s - 1);
b = zeros(2, m, s);
d = zeros(2, s);
c = zeros(s - 1, 1);

% -- A Matrix ------------------------------------------------------------------------

% -- b Vector ------------------------------------------------------------------------
% \varphi_1(hL) * F(Y1)
b(:, 2, 1) = [1.0,  1.0];

% -- d Vector ------------------------------------------------------------------------
d(:, 1) = [1.0,  1.0]; % exp(h L) y_n

% -- c Vector ------------------------------------------------------------------------

end

