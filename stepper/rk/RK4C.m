function [A, b, c] = RK4C()
% Classical RK4 integrator
 
A = [
    1 / 2,              0,              0;
    0,                  1/2             0;
    0,                  0,              1
];
b = [ 1 / 6,    1 / 3,    1 / 3,    1 / 6 ];
c = [1/2 1/2 1];

end