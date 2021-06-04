function [A_im, b_im, A_ex, b_ex, c] = IMRK1()
%IMRK1 Forward-backward Euler

% -- implicit coefficients -----------------------------------------------------------
A_im = [0 0; 0 1];
b_im = [0 1];

% -- explicit coefficients -----------------------------------------------------------
A_ex = [0 0; 1 0];
b_ex = [1 0];

% -- c Vector ------------------------------------------------------------------------
c    = [0 1];

end

