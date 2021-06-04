function [R] = rIMRK4(z1, z2_vec)
%rERK4 stability function for IMEX ARK4(3)6L[2]SA from "C. A. Kennedy and M. H. Carpenter, "Additive Runge-Kutta 
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
%INITETDRK4 Initializes coefficients for ERK1 for scaler Lambda

A_im = zeros(6);
A_im(2,1) =   0.25;
A_im(2,2) =   0.25;
A_im(3,1) =   8611.0        / 62500.0;
A_im(3,2) = - 1743.0        / 31250.0;
A_im(3,3) =   0.25;
A_im(4,1) =   5012029.0     / 34652500.0;
A_im(4,2) = - 654441.0      / 2922500.0;
A_im(4,3) =   174375.0      / 388108.0;
A_im(4,4) =   0.25;
A_im(5,1) =   15267082809.0 / 155376265600.0;
A_im(5,2) = - 71443401.0    / 120774400.0;
A_im(5,3) =   730878875.0   / 902184768.0;
A_im(5,4) =   2285395.0     / 8070912.0;
A_im(5,5) =   0.25;     
A_im(6,1) =   82889.0       / 524892.0;
A_im(6,2) =   0.0;
A_im(6,3) =   15625.0       / 83664.0;
A_im(6,4) =   69875.0       / 102672.0;
A_im(6,5) = - 2260.0        / 8211.0;
A_im(6,6) =   0.25;
      
b_im = [82889.0 / 524892.0, 0.0, 15625.0 /  83664.0, 69875.0 / 102672.0, - 2260.0 / 8211.0, 0.25];  

A_ex = zeros(6);
A_ex(2,1) =   0.5;
A_ex(3,1) =   13861.0          / 62500.0;
A_ex(3,2) =   6889.0           / 62500.0;
A_ex(4,1) = - 116923316275.0   / 2393684061468.0;
A_ex(4,2) = - 2731218467317.0  / 15368042101831.0;
A_ex(4,3) =   9408046702089.0  / 11113171139209.0;
A_ex(5,1) = - 451086348788.0   / 2902428689909.0;
A_ex(5,2) = - 2682348792572.0  / 7519795681897.0;
A_ex(5,3) =   12662868775082.0 / 11960479115383.0;
A_ex(5,4) =   3355817975965.0  / 11060851509271.0;
A_ex(6,1) =   647845179188.0   / 3216320057751.0;
A_ex(6,2) =   73281519250.0    / 8382639484533.0;
A_ex(6,3) =   552539513391.0   / 3454668386233.0;
A_ex(6,4) =   3354512671639.0  / 8306763924573.0;
A_ex(6,5) =   4040.0           / 17871.0;

b_ex = b_im;

end