

function fit_PD(MaxTime,n_epochs, make_plot)

if nargin < 3
	make_plot = false;
end


load('all_PD')
V0 = all_PD;
dV0 = NaN*V0;
T = 10;
for i = 1:size(all_PD,2)
	V0(:,i) = filtfilt(ones(T,1),T,V0(:,i));
	dV0(:,i) = [NaN; diff(V0(:,i))];
end

% prepare data to be sent to cost function
data = struct;
for i = 1:size(all_PD,2)

	[spike_amplitudes, spike_peaks, minima_bw_spikes] = measurePDmetrics(V0(:,i));

	spike_amplitude_range = [nanmin(spike_amplitudes) nanmax(spike_amplitudes)];


	smallest_spike = min(spike_peaks) - .5;


	metrics =  xtools.V2metrics(V0(:,i),'spike_threshold',smallest_spike,'sampling_rate',10);

	data(i).max_V = max(V0(:,i));
	data(i).min_V = min(V0(:,i));
	data(i).burst_period = metrics.burst_period;
	data(i).duty_cycle = metrics.duty_cycle_mean;
	data(i).V0 = V0(:,i);
	data(i).dV0 = dV0(:,i);

	data(i).spike_amplitude_range = spike_amplitude_range;

	data(i).V_bw_spikes_range = [min(minima_bw_spikes) max(minima_bw_spikes)];
	data(i).spike_peaks = [min(spike_peaks) max(spike_peaks)];

	data(i).n_spikes_per_burst = metrics.n_spikes_per_burst_mean;

end


% set up the neuron
x = xolotl.examples.TwoCompartmentSTG;

% set up the fitter
p = xfit;
p.x = x;

p.parameter_names = [x.find('*gbar');'Axon.len'; 'CellBody.len'; 'CellBody.radius'; 'CellBody.CalciumMech.f'; 'CellBody.CalciumMech.tau_Ca'];
%      Axon                Soma  
%      A   Kd    L    NaV   A   CaS  CaT  H    KCa  Kd   L
p.lb = [10,  100, 0,  100, 10,  1,   1,   1,   10,  10,  0,  .1, .01, .01, 1,  10];
p.ub = [1e3, 2e3, 10, 3e3, 1e3, 100, 100,  10, 1e3, 2e3, 10,  5,  1,  1,  20, 500];




p.data = data;

% debug
p.options.UseParallel = true;

p.sim_func = @metricsCost;



% run forever, churning out neurons
N = 1e4;
all_params = NaN(N,length(p.ub));
all_cost = NaN(N,1);

file_name = [corelib.getComputerName '_PD_models.mat'];

if exist(file_name)
	load(file_name)
	start_idx = find(~isnan(all_cost),1,'last')+1;
else
	start_idx = 1;
end

p.options.MaxTime = MaxTime;
p.options.Display = 'iter';


for i = start_idx:N
	disp(['Starting with random seed #' strlib.oval(i)])
	try

		% configure the data
		p.data = data(rem(i,length(data))+1);

		p.seed = rand(length(p.seed),1).*(p.ub - p.lb) + p.lb;
		for j = 1:n_epochs
			p.fit;
		end

		% save
		X = p.seed;
		x.set(p.parameter_names,X);
		this_cost = p.sim_func(x, p.data);

		if isnan(this_cost)
			continue
		end

		if make_plot
			figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
			x.integrate;
			V = x.integrate;
			time = (1:length(V))*x.dt*1e-3;
			plot(time,V(:,2),'k')
			set(gca,'XLim',[0 3])


		end

		if this_cost < 4


			disp(['Final Cost for this seed is ' strlib.oval(this_cost)])

			all_params(i,:) = X;
			all_cost(i) = this_cost;

			save(file_name,'all_params','all_cost')


		else
			disp('Cost too high, skipping...')
		end

	catch
		disp('Something went wrong. Ouch. ')
	end


end
