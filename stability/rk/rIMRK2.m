function [R] = rIMRK2(z1, z2_vec)
%rIMRK2 stability function for IMEX (2,3,2) methods from U. M. Ascher, S. J. Ruuth, R. J. Spiteri, "Implicit-Explicit
% Runge-Kutta methods for time-dependent partial differential equations" 
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
% RETURNS
%   R    - value of stability function

[A_im, b_im, A_ex, b_ex] = initIMEXCoefficients();
R = rIMRK(z1, z2_vec, A_im, b_im, A_ex, b_ex);
end

function [A_im, b_im, A_ex, b_ex] = initIMEXCoefficients()
%INITIMEXCOEFFICIENTS initialized butcher tablaues for IMEX method

gamma  = (2 - sqrt(2))/2;
delta  = -2*sqrt(2)/3;

A_im = zeros(3);
A_im(2,2) = gamma;
A_im(3,2) = 1 - gamma;
A_im(3,3) = gamma;
b_im = [0, 1 - gamma, gamma];

A_ex = zeros(3);
A_ex(2,1) = gamma;
A_ex(3,1) = delta;
A_ex(3,2) = 1 - delta;
b_ex = [0, 1 - gamma, gamma];

end