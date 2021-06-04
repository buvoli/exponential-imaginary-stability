function [data_raw] = TDStabilityRegionData(amp, z1_r, z2_r, z1_angle, z2_angle)
%TDSTABILITYDATA produces the raw data for an two-dimensional stability
%plot where z1 and z2 are complex numbers with fixed angle and variable
%radius.
%   amp   (handle)  - function of two arguments @(z1 - scalar, z2 - vector) producing amp factors
%   z1_r  (vector)  - radius of z1 component. 
%   z2_r  (vector)  - radius of z2 component. 
%   z1_angle (real) - angle of z1 component
%   z2_angle (real) - angle of z2 component

num_z1 = length(z1_r);
num_z2 = length(z2_r);
data_raw = zeros(num_z2, num_z1);
for i = 1 : num_z1
    data_raw(:,i) = amp(z1_angle * z1_r(i), z2_angle * z2_r);
end
end