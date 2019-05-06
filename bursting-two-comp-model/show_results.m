% this script shows the best solutions we found

x = makeSTGNeuron;


x.dt = .1;
x.t_end = 3e3;


parameter_names = [x.find('*gbar');'Axon.len'; 'CellBody.len'; 'CellBody.radius'; 'CellBody.CalciumMech.f'; 'CellBody.CalciumMech.tau_Ca'];

% load the result
load('PD_models.mat')

N = find(~isnan(all_cost),1,'last');
show_these = veclib.shuffle(1:N);

show_these = N:-1:1;

figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
ph = plot(NaN,NaN,'k');

for i = 1:length(show_these)
	x.set(parameter_names,all_params(show_these(i),:))

	x.reset;
	x.integrate;
	x.integrate;
	V = x.integrate;
	time = (1:length(V))*x.dt*1e-3;
	ph.XData = time;
	ph.YData = V(:,2);

	title(show_these(i))

	pause(1)

end