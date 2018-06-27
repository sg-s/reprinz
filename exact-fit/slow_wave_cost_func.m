function C = slow_wave_cost_func(x,data)

inf_cost = 1e3;
C = inf_cost;


x.CellBody.V = -54;

V = x.integrate;


C = sum((V - data.slow_wave).^2);