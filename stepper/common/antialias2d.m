function [yh] = antialias2d(yh)
%ANTIALIAS2D Summary of this function goes here
% performs 2d antialias 2/3 filter on vector yh
% zeros out modes (Nx/2 - b) ... (Nx/2) and (-Nx/2 + 1) ... (-Nx / 2 + b)
% === Parameters ===
% yh - 1d solution in fourier space. Wave
% === Output ===
% yh - anti aliased 1d solution in fourier space
    
    NpX = size(yh, 1);
    NpY = size(yh, 2);
    
    aX  = floor(NpX / 2.0);
    bX  = floor(NpX / 6.0);

    aY  = floor(NpY / 2.0);
    bY  = floor(NpY / 6.0);
      
    yh( aX + 1 - bX : aX + 1 + bX, : ) = 0.0;
    yh( : , aY + 1 - bY : aY + 1 + bY ) = 0.0;
end