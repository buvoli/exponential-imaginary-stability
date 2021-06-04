function [h, e] = HTrap(fh, params)
%H Summary of this function goes here
%   Detailed explanation goes here

Lv = params.Lv; Nv = params.Nv;
vs = linspace(-Lv / 2, Lv / 2, Nv + 1); vs(end) = [];

f = ifft(reshape(fh, [params.Nv params.Nx]), [], 2);

dv = params.Lv / params.Nv;
dx = params.Lx / params.Nx;

I1 = dx * dv * sum(sum(bsxfun(@times, vs(:).^2, f)));
I2 = dx * sum(abs(params.E(fh, params).^2));
h  = I1 / 2 + I2 / 2;
e  = I2;
end