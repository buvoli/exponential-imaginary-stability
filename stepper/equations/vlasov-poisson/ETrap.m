function E = ETrap(fh, params)
%E Computes the energy E which satisfies
%  E_x = \int_{R} f(t,x,v) dx - 1
% PARAMETERS
%  fh - soution vector where x is in fourier space and k is in physical
%  space.
% Returns
%  E - energy in physical space
dv = params.Lv / params.Nv;
I = dv * sum(reshape(fh, [params.Nv params.Nx]), 1);
I(:,1) = I(:,1) - params.Nv; % add fft(-1 * ones(size(I)))
E = ifft(params.IDX(1,:) .* I);
end