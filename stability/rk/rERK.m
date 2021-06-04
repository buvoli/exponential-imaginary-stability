function [R] = rERK(z2_vec, d, A, b)
% rERK computes stability function for the exponential runge kutta method
%
%      Y_i = d(i) y_n + \sum_{j=1}^{i-1} A(i,j) F(c_j, Y_j)    i = 1, ..., s
%      y_{n+1} = d(i) y_n + \sum_{j=1}^s b(i) F(c_j, Y_j)
%
% PARAMETERS
%   d     (vector (s+1)x1) - coefficient for y_n 
%   A     (matrix sxs) - coefficients for combining stage derivatives
%   b     (vector sx1) - output coefficient for combining stage derivatives
%   z2    (scalar or vector) - explicit term z_2 = h * \lambda_2
% RETURNS
%   R    - value of stability function

sz = size(z2_vec);
z2_vec = reshape(z2_vec, 1, length(z2_vec)); % ensure row vector

s = size(A, 1);  % num stages
Y = zeros(s, length(z2_vec)); % stage values

for i = 1 : s
    Y(i,:) = d(i);
    for j = 1 : (i-1)
        Y(i,:) = Y(i,:) + A(i,j) * z2_vec .* Y(j,:);
    end
end

R = d(s+1) * ones(1, length(z2_vec));
for i = 1 : s
    R = R + (b(i) * z2_vec .* Y(i,:));
end
R = reshape(R, sz); % resize to same dimensions as z2_vec
end