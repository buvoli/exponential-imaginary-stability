function [fh, data_raw] = TDRStabilityRegionPlot(amp, z1_r, z2_r, z1_angle, z2_angle, rho, options)
%TDStabilityRegionPlot produces a two dimensional repartitioned (TDR) plot
%   amp    (handle)  - function of two arguments @(z1 - scalar, z2 - vector) producing amp factors
%   z1_r   (vector)  - radius of z1 component. 
%   z2_r   (vector)  - radius of z2 component. 
%   z1_angle  (real) - angle of z1 component
%   z2_angle  (real) - angle of z2 component
%   rho (real) - desired angle delta

if(nargin == 3)
    options = struct();
end

data_raw = TDRStabilityRegionData(amp, z1_r, z2_r, z1_angle, z2_angle, rho);
fh       = TDStabilityRegionPlotter(data_raw, z1_r, z2_r, options);
end