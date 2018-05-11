
x = make_bursting_soma();

parameter_names = {'Soma.ACurrent.gbar', 'Soma.CaS.gbar', 'Soma.CaT.gbar', 'Soma.HCurrent.gbar', 'Soma.KCa.gbar' , 'Soma.Kd.gbar' , 'Soma.Leak.gbar' };

[~,show_these] = sort(all_cost);


for i = 1:length(show_these)


	
	x.set(parameter_names,all_g(:,show_these(i)));

	bursting_soma_cost_func(x);
	title(oval(show_these(i)))
	drawnow;
	pause(2)
	close all
	clc
end


