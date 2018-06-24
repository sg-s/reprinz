function C = exact_fit_cost_func(x)


% copy params
g = x.get('Soma*gbar');
x.set('Neurite*gbar',g);
x.Soma.NaV.gbar = 0;

x.Soma.V = -54;
x.Neurite.V = -54;
I_clamp  = x.integrate();

C = mean(abs(I_clamp(20e3:end)));

