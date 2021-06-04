init;
Nt = 800; [ts, ys, tcpu] = ERK(LF(pars), NF, tspan, y0, Nt, struct('coeffGenerator', @ERK4, 'parameters', pars));
figure(); surf(xs, vs, fliplr(to_physical(ys(:,end)))); view([90 90]); axis tight; shading interp; colorbar; caxis([0 0.4]);

%% Energy Drift Test
ns = size(ys,2);
en = zeros(ns, 2);

for i = 1 : ns
    [en(i,1), en(i,2)] = HTrap(ys(:,i), pars);
end
en(:,1) = (en(:,1) - en(1,1)) / en(1,1);

figure(); semilogy(ts, abs(en(:,1))); title('Relative energy Drift');
figure(); plot(ts, en(:, 2)); title('energy \|E\|^2');