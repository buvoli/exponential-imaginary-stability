function [t0, y0, ts, ys] = rk4(f, tspan, y0, Nt, options)
%RK4 Classical RK4 integrator
%   Detailed explanation goes here
 
t0 = tspan(1);
h  = diff(tspan) / Nt;
pars = options.parameters;

if(nargout >= 3)
    store_data = true;
    if(isfield(options, 'max_ts_to_store'))
        max_ts_to_store = options.max_ts_to_store;
    else
        max_ts_to_store = max(2,min(5000000/length(y0),1000));
    end
    skip_rate = ceil(Nt/max_ts_to_store);
    ys = zeros(length(y0),ceil(Nt/skip_rate)+1);
    ts = zeros(size(ys,2),1);
    ys(:,1) = y0; save_count = 2;
    ts(1) = tspan(1); t0 = tspan(1);
else
    store_data = false;
end
 
for i = 1:Nt
    
    k1 = f(t0, y0, pars)*h;
    k2 = f(t0 + h/2, y0 + k1/2, pars) * h;
    k3 = f(t0 + h/2, y0 + k2/2, pars) * h;
    k4 = f(t0 + h, y0 + k3, pars) * h;
    y0 = y0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    t0 = t0 + h;
    
    % Save Data
    if(store_data && (mod(i,skip_rate) == 0 || i==Nt))
        ys(:,save_count) = y0;
        ts(save_count,:) = t0;
        save_count = save_count + 1;
    end 
        
end
 
end