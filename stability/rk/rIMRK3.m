function [R] = rIMRK3(z1, z2_vec)
%rIMRK1 stability function for IMEX ARK3(2)4L[2]SA from "C. A. Kennedy and M. H. Carpenter, "Additive Runge-Kutta 
% schemes for convection-diffusion-reaction equations", (2003). 
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

A_ex = zeros(4);
A_im = zeros(4);

A_ex(2,1) =   1767732205903.0  / 2027836641118.0;
A_ex(3,1) =   5535828885825.0  / 10492691773637.0;
A_ex(3,2) =   788022342437.0   / 10882634858940.0;
A_ex(4,1) =   6485989280629.0  / 16251701735622.0;
A_ex(4,2) = - 4246266847089.0  / 9704473918619.0;
A_ex(4,3) =   10755448449292.0 / 10357097424841.0;

b_ex = [
    1471266399579.0  / 7840856788654.0
    - 4482444167858.0 / 7529755066697.0
    11266239266428.0 / 11593286722821.0
    1767732205903.0 / 4055673282236.0
];

A_im(2,1) =   1767732205903.0  / 4055673282236.0;
A_im(2,2) =   1767732205903.0  / 4055673282236.0;
A_im(3,1) =   2746238789719.0  / 10658868560708.0;
A_im(3,2) = - 640167445237.0   / 6845629431997.0;
A_im(3,3) =   1767732205903.0  / 4055673282236.0;
A_im(4,1) =   1471266399579.0  / 7840856788654.0;
A_im(4,2) = - 4482444167858.0  / 7529755066697.0;
A_im(4,3) =   11266239266428.0 / 11593286722821.0;
A_im(4,4) =   1767732205903.0  / 4055673282236.0;

b_im = b_ex;

end