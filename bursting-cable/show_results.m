
x = make_realistic_neuron();
x.sha1hash;

parameter_names = {'Soma.ACurrent.gbar', 'Soma.CaS.gbar', 'Soma.CaT.gbar', 'Soma.HCurrent.gbar', 'Soma.KCa.gbar' , 'Soma.Kd.gbar' , 'Soma.NaV.gbar' };

show_these = shuffle(find(~isnan(all_cost)));

for i = 1:length(show_these)


	
	x.set(parameter_names,all_g(:,show_these(i)));

	realistic_cost_func(x);
	title(oval(show_these(i)))
	drawnow;
	pause(2)
	close all
	clc
end


