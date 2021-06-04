function ns = N(~,fhs,params)
%N Returns ODE nonlinear function
% fh - fourier transform of f (in x only) -- f(k,v,t)

num_fh = size(fhs,2);
if( num_fh == 0) % runs faster without loop for single inputs
    ns = Nv(fhs, params);
    return;
else
    ns = zeros(params.Nv * params.Nx, num_fh);
    for i = 1 : num_fh
        ns(:,i) = Nv(fhs(:,i), params);
    end
end

end

function n = Nv(fh,params)

    fh = reshape(fh,[params.Nv, params.Nx]);
    E  = E_trap(fh, params);

    switch(params.vdes)
        case 'CD'
            f_v = params.DV * ifft(fh, [], 2);
        case 'FS'
           f_v = ifft(params.DV .* fft(ifft(fh, [], 2), [], 1), [], 1);
        case 'W5'
           f_v = DV_Weno5(ifft(fh, [], 2), E, params);
    end

    
    n = fft(-1 * E .* f_v, [], 2) + repartitionTerm(fh,params);
    
    if(params.antialias)
        n = ifft(antialias2d(fft(n,[],1)),[],1);
    end
    
    n = reshape(n, [params.Nv * params.Nx, 1]);

end

function E = E_trap(fh, params)
% fh - fourier transform of f (in x only) -- f(k,v,t)
% returns E in physical x space

% -- naive implementation -------------------------------------------------
% f = ifft(fh, [], 2);
% dx = params.Lv / params.Nv;
% I = -1 + dx * sum(f, 1);
% E = ifft(params.IDX .* fft(I,[],2), [], 2);

% -- faster equivalent implementation -------------------------------------
dv = params.Lv / params.Nv;
I = dv * sum(fh, 1);
I(:,1) = I(:,1) - params.Nv; % add fft(-1 * ones(size(I)))
E = ifft(params.IDX .* I, [], 2);

end

% == WENO Functions =======================================================

function f_v = DV_Weno5(f, E, params)

% evaluate smoothness indicator stencils
SI1 = params.DV{1,1} * f;
SI2 = params.DV{1,2} * f;
SI3 = params.DV{1,3} * f;
SI4 = params.DV{1,4} * f;

% evaluate derivative stencils
DS1 = params.DV{2,1} * f;
DS2 = params.DV{2,2} * f;
DS3 = params.DV{2,3} * f;
DS4 = params.DV{2,4} * f;

fp_v = upwindDv(SI1, SI2, SI3, SI4, DS1, DS2, DS3, DS4, params);
fm_v = downwindDv(SI1, SI2, SI3, SI4, DS1, DS2, DS3, DS4, params);
f_v  = (real(E) > 0) .* fp_v + (real(E) < 0) .* fm_v;
return;

end

function fp_v = upwindDv(SI1, SI2, SI3, SI4, DS1, DS2, DS3, DS4, params) % upwind WENO

    Nv = params.Nv;

    B0 = B( 13 / 12, SI1, -1, (1/4), SI2, -1, Nv);
    B1 = B( 13 / 12, SI1,  0, (1/4), SI4,  0, Nv);
    B2 = B( 13 / 12, SI1,  1, (1/4), SI3,  1, Nv);

    [w0, w1, w2] = W(B0,B1,B2);

    fp_jp = F(w0, w1, w2, DS1, -1, DS2, 0, DS3, 1, Nv);
    fp_jm = fp_jp(pInds(Nv, -1), :);

    dv = params.Lv / params.Nv;
    fp_v = (fp_jp - fp_jm) / dv;

end

function fm_v = downwindDv(SI1, SI2, SI3, SI4, DS1, DS2, DS3, DS4, params) % downwind WENO

    Nv = params.Nv;

    B0 = B( 13 / 12, SI1, 2, (1/4), SI3,  2, Nv);
    B1 = B( 13 / 12, SI1, 1, (1/4), SI4,  1, Nv);
    B2 = B( 13 / 12, SI1, 0, (1/4), SI2,  0, Nv);

    [w0, w1, w2] = W(B0,B1,B2);

    fm_jp = F(w0, w1, w2, DS2, 0, DS3, 1, DS4, 2, Nv);
    fm_jm = fm_jp(pInds(Nv, -1), :);
    
    dv = params.Lv / params.Nv;
    fm_v = (fm_jp - fm_jm) / dv;

end

function b = B(c1, S1, shift1, c2, S2, shift2, N)   
    ind1 = pInds(N, shift1);
    ind2 = pInds(N, shift2);
    b = c1 * S1(ind1, :) .^2 + c2 * S2(ind2, :) .^2;
end

function [w0,w1,w2] = W(B0, B1, B2)

    % parameters
    gamma0 = 1 / 10;
    gamma1 = 6 / 10;
    gamma2 = 3 / 10;
    epsilon = 1e-6;

    A0 = gamma0 ./ (epsilon + B0) .^ 2;
    A1 = gamma1 ./ (epsilon + B1) .^ 2;
    A2 = gamma2 ./ (epsilon + B2) .^ 2;
    Atot = (A0 + A1 + A2);
    
    w0 = A0 ./ Atot;
    w1 = A1 ./ Atot;
    w2 = A2 ./ Atot;
    
end

function f = F(W0, W1, W2, S1, shift1, S2, shift2, S3, shift3, N)   
    ind1 = pInds(N, shift1);
    ind2 = pInds(N, shift2);
    ind3 = pInds(N, shift3);
    f = W0 .* S1(ind1, :) + W1 .* S2(ind2, :) + W2 .* S3(ind3, :);
end

function i = pInds(N, shift) % periodic indices
    i = mod((1:N) + shift - 1, N) + 1;
end

% == Repartitioning Functions =============================================

function D = repartitionTerm(fh, params, V)

    D = 0;
    
    flag_trans   = params.epsilon ~= 0;
    flag_rot     = params.rho ~= 0;
    flag_overrot = params.rhouxx ~= 0;
    
    if( flag_trans || flag_rot || flag_overrot )
        xs = linspace(0, params.Lx, params.Nx + 1); xs(end) = [];
        vs = linspace(-params.Lv / 2, params.Lv / 2, params.Nv + 1); vs(end) = [];
        [~, V] = meshgrid(xs,vs);
    end
    
    if( flag_trans ) % translation
        D = D + abs(V) .* params.epsilon .* fh;
    end

    if( flag_rot ) % rotation
        D = D - abs(V) .* abs(params.DX) / tan(pi/2 + params.rho) .* fh;
    end
    
    if( flag_overrot ) % over-rotation
        D = D + abs(V) .* params.DXX / tan(pi/2 + params.rhouxx) .* fh;
    end

end