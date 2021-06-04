function cs = L(params)
%L Returns coefficients for linear ODEs

xs = linspace(0, params.Lx, params.Nx + 1); xs(end) = [];
vs = linspace(-params.Lv / 2, params.Lv / 2, params.Nv + 1); vs(end) = [];
[~, V] = meshgrid(xs,vs); % V = repmat(vs(:), pars.Nx)

% Linear Operator
L = -V .* params.DX + repartitionMatrix(params, V);

if(params.antialias)
    L = ifft(antialias2d(fft(L,[],1)),[],1);
end

% Reshape to row vector format
cs = reshape(L, params.Nv*params.Nx, 1);
end

function D = repartitionMatrix(params, V)

    D = 0;
   
    if(params.epsilon ~= 0) % translation
        D = D - abs(V) * params.epsilon;
    end

    if(params.rho ~= 0) % rotation
        D = D + abs(V) .* abs(params.DX) / tan(pi/2 + params.rho);
    end
    
    if(params.rhouxx ~= 0) % over-rotation
        D = D - abs(V) .* params.DXX / tan(pi/2 + params.rhouxx);
    end

end