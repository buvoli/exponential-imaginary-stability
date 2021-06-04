function [R] = rIMRK(z1, z2_vec, A_im, b_im, A_ex, b_ex)
% rRK computes stability function for the IMEX runge kutta method
%
%      Y_i = y_n + \sum_{j=1}^{i-1} A(i,j) F(c_j, Y_j)    i = 1, ..., s
%      y_{n+1} = d(i) y_n + \sum_{j=1}^s b(i) F(c_j, Y_j)
%
% PARAMETERS
%   d      (vector (s+1)x1) - coefficient for y_n 
%   A      (matrix sxs) - coefficients for combining stage derivatives
%   b      (vector sx1) - output coefficient for combining stage derivatives
%   z1     (scalar) - implicit term z_1 = h * \lambda_1
%   z2_vec (scalar or vector) - explicit term z_2 = h * \lambda_2
% RETURNS
%   R    - stability function

sz = size(z2_vec);
z2_vec = reshape(z2_vec, 1, length(z2_vec)); % ensure row vector

s = size(A_im, 1);  % num stages
Y = zeros(s, length(z2_vec)); % stage values

for i = 1 : s
    Y(i,:) = 1;
    for j = 1 : (i-1)
        Y(i,:) = Y(i,:) + A_im(i,j) * z1 .* Y(j,:) + A_ex(i,j) * z2_vec .* Y(j,:);
    end
    Y(i,:) = Y(i,:) / (1 - z1 * A_im(i,i));
end

R = ones(1, length(z2_vec));
for i = 1 : s
    R = R + (b_im(i) * z1 + b_ex(i) * z2_vec) .* Y(i,:);
end
R = reshape(R, sz); % resize to same dimensions as z2_vec
end