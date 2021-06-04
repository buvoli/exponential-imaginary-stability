function [DX, DV, IDX, DXX] = DM(params)
%DM Returns differentation Matrices & Inverse Lapalacian for 2D periodic Domain

%set up wavenumbers
wn_x = wavenumbers(params.Nx, params.Lx);

% set up Linear L Matrix
m_x = spdiags(wn_x,0,params.Nx,params.Nx);
O   = ones(params.Nv,params.Nx);

% Build Differentiation coefficients
DX  = O*(1i*m_x);
DXX  = O*(1i*m_x).^2; % for repartitioning

switch(params.vdes)
    case 'CD'
        DV = DV_CD(params);
    case 'FS'
        DV = DV_F(params);
    case 'W5'
        DV = DV_WENO5(params);
    otherwise
        DV = DV_CD(params);
end

% Build Inverse Integral coefficients 
IDX = 1 ./ DX;
IDX(:,1) = 0; %Set (0,0) mode to zero

end

function wn = wavenumbers(N, L)
    if(mod(N,2) == 0)
        wn = [0:N/2 -N/2+1:-1].' * (2*pi/(L));
    else
        wn = [0:(N-1)/2, -(N-1)/2:-1].' * (2*pi/L);
    end
end

function DV = DV_CD(params)
% 2nd-order centered difference for u_v
    dv = params.Lv / params.Nv; 
    w  = [-1 0 1] / (2 * dv);
    DV = sTriCirc(w, params.Nv);
end

function DV = DV_F(params)
% Fourier spectral derivative operator for u_v
    
    wn_v = wavenumbers(params.Nv, params.Lv);%[0:params.Nv/2 -params.Nv/2+1:-1].' * (2*pi/(params.Lv));
    m_v = spdiags(wn_v,0,params.Nv,params.Nv);
    O   = ones(params.Nv, params.Nx);
    DV  = (1i*m_v)*O;
end

function DV = DV_WENO5(params)
% Weno5 spectral derivative operator for u_v
    
    % Smoothness Indicator Stencils
    SS1 = sTriCirc([ 1, -2,  1], params.Nv);
    SS2 = sTriCirc([ 1, -4,  3], params.Nv);
    SS3 = sTriCirc([ 3,  4,  1], params.Nv);
    SS4 = sTriCirc([ 1,  0, -1], params.Nv);
    
    % Derivative Stencils
    DS1 = sTriCirc([ 2,  -7,  11] / 6, params.Nv);
    DS2 = sTriCirc([-1,   5,  2] / 6, params.Nv);
    DS3 = sTriCirc([ 2,   5,  -1] / 6, params.Nv);
    DS4 = sTriCirc([ 11, -7,  2] / 6, params.Nv);
    
    DV = {
        SS1, SS2, SS3, SS4;
        DS1, DS2, DS3, DS4
    };

end

function TC=sTriCirc(w, N)
% creates a sparse, tridagonal circulant matrix

e   = ones(N, 1);
TC  = spdiags([w(1) * e, w(2) * e, w(3) * e], -1:1, N, N);
TC(1, N) = w(1);
TC(N, 1) = w(3);
end