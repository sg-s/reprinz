
x = make_realistic_neuron();
x.sha1hash;
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = {'Soma.ACurrent.gbar', 'Soma.CaS.gbar', 'Soma.CaT.gbar', 'Soma.HCurrent.gbar', 'Soma.KCa.gbar' , 'Soma.Kd.gbar' , 'Soma.NaV.gbar' };

show_these = (find(~isnan(all_cost)));

for i = 1:length(show_these)


	
	x.set(p.parameter_names,all_g(:,show_these(i)));

	realistic_cost_func(x);
	title(oval(i))
	drawnow;
	pause(2)
	close all
end


