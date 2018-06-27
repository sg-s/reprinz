
x = make_bursting_soma();

parameter_names = {'CellBody.ACurrent.gbar', 'CellBody.CaS.gbar', 'CellBody.CaT.gbar', 'CellBody.HCurrent.gbar', 'CellBody.KCa.gbar' , 'CellBody.Kd.gbar' , 'CellBody.Leak.gbar' };

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


