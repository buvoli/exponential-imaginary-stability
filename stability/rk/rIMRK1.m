function [R] = rIMRK1(z1, z2_vec)
%rIMRK1 stability function for IMEX forward/backward euler pair
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

A_im = [0 0; 0 1];
b_im = [0 1];
A_ex = [0 0; 1 0];
b_ex = [1 0];
end