function C = exact_fit_cost_func(x, I_ext, V_clamp)


x.AB.V = V_clamp(1);
I_clamp  = x.integrate([],V_clamp');

C = mean(abs(I_clamp(20e3:end)));
