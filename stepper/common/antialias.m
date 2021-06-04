function [yh] = antialias(yh, dim)
%ANTIALIAS Summary of this function goes here
% performs 1d antialias 2/3 filter on vector yh
% zeros out modes (Nx/2 - b) ... (Nx/2) and (-Nx/2 + 1) ... (-Nx / 2 + b)
% === Parameters ===
% yh - 1d solution in fourier space. Wave
% === Output ===
% yh - anti aliased 1d solution in fourier space
    
    if(isvector(yh))
        Np = length(yh);
        a  = floor(Np / 2.0);
        b  = floor(Np / 6.0);
        yh( a + 1 - b : a + 1 + b) = 0.0;
    else
        yh_size = size(yh);        
        Np = yh_size(dim);
        a  = floor(Np / 2.0);
        b  = floor(Np / 6.0);
        inds = arrayfun(@(s) 1:s, yh_size, 'UniformOutput', false);
        inds{dim} = a + 1 - b : a + 1 + b;
        yh(inds{:}) = 0.0;        
    end
end

