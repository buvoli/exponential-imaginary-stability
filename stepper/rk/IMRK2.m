function [A_im, b_im, A_ex, b_ex, c] = IMRK2()
%IMRK2 Ascher-Ruuth-Spiteri

gamma  = (2 - sqrt(2))/2;
delta  = -2*sqrt(2)/3;

% -- implicit coefficients -----------------------------------------------------------
A_im = zeros(3);
A_im(2,2) = gamma;
A_im(3,2) = 1 - gamma;
A_im(3,3) = gamma;
b_im = [0, 1 - gamma, gamma];

% -- explicit coefficients -----------------------------------------------------------
A_ex = zeros(3);
A_ex(2,1) = gamma;
A_ex(3,1) = delta;
A_ex(3,2) = 1 - delta;
b_ex = [0, 1 - gamma, gamma];

% -- c Vector ------------------------------------------------------------------------
c    = [0 gamma 1];

end