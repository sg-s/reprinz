
clear x
x = make_realistic_cable();
x.md5hash;

parameter_names = {'CellBody.ACurrent.gbar', 'CellBody.CaS.gbar', 'CellBody.CaT.gbar', 'CellBody.HCurrent.gbar', 'CellBody.KCa.gbar' , 'CellBody.Kd.gbar' , 'CellBody.NaV.gbar' };

show_these = shuffle(find(~isnan(all_cost)));

show_these = 615;

for i = 1:length(show_these)


	
	x.set(parameter_names,all_g(:,show_these(i)));

	realistic_cost_func(x);
	title(oval(show_these(i)))
	drawnow;
	pause(2)
	close all
	clc
end


