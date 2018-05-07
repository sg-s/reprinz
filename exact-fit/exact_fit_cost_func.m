function C = exact_fit_cost_func(x, I_ext, V_clamp)


% copy params
g = x.get(x.find('Soma*gbar'));
x.set(x.find('Neurite*gbar'),g);
x.Soma.NaV.gbar = 0;

x.Soma.V = V_clamp(1);
x.Neurite.V = V_clamp(1);
% I_clamp  = x.integrate([],V_clamp');

% C = mean(abs(I_clamp(20e3:end)));

V = x.integrate;
C = sum((V(:,2) - V_clamp).^2);