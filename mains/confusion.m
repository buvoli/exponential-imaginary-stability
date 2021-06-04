% cannot explain why following is stable (set rho to pi/8)

init; [ts,ys] = EPBM_PMFCmS(LF(pars), NF, 2 * tspan, y0, 2 * 800, struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 2, 'parameters', pars)); mesh(sfilter(ys))
init; [ts,ys] = EPBM_PMFCmS(LF(pars), NF, 2 * tspan, y0, 2 * 800, struct('z', lobpts(4), 'b', -1 * ones(4), 'mS', [], 'kappa', 1, 'parameters', pars)); mesh(sfilter(ys))

% using rho = pi/64, and tspan = [0 80]
init; [ts,ys] = EPBM_PMFCmS(LF(pars), NF, 2 * tspan, y0, 4 * 800, struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 0, 'parameters', pars)); mesh(sfilter(ys))
% however, it gets worse if kappa > 0?



%%


z1s = linspace(0, 20,200);
z2s = linspace(-2,2,100);
a = zeros(length(z1s), length(z2s));
a_new = zeros(length(z1s), length(z2s));
rho0 = pi / 8;

for i = 1 : length(z1s)
    shift = z1s(i)^(3/3) / tan(pi/2 + rho0);
    a(i,:) = ampETDPBM(1i * z1s(i) + shift, 1i * z2s(:) - shift, struct('z', lobpts(4), 'b', -1 * ones(4), 'mS', [], 'kappa', 3, 'alpha', 2));%struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 2, 'alpha', 2)); 
    %a(i,:) = ampETDPBM(1i * z1s(i) + shift, 1i * z2s(:) - shift, struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 1, 'alpha', 1)); 
%    a(i,:) = abs(rERK4(1i * z1s(i) + shift, 1i * z2s(j) - shift));
%     for j = 1 :length(z2s)
%         shift = z1s(i) / tan(pi/2 + rho0);
%         a(i,j) = ampETDPBM(1i * z1s(i) + shift, 1i * z2s(j) - shift, struct('z', [-1 ; legpts(4)], 'b', -1 * ones(5), 'mS', 1, 'kappa', 0, 'alpha', 2));
%         %a(i,j) = abs(rERK4(1i * z1s(i) + shift, 1i * z2s(j) - shift));
%     end
end

figure(100);

%a(a > 1) = NaN;
surf(z2s, z1s, a); shading interp; view([90 90]); hold on;
[~, cnt] = contour(z2s, z1s, a, [0 1], 'k', 'LineWidth', 2); hold off;
cnt.ZLocation = 1;

