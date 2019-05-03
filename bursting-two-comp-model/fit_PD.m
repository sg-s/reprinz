

function fit_PD(MaxTime,n_epochs)

load('../exact-fit/isolated_PD')
V0 = V;
T = 1e3;
V = filtfilt(ones(T,1),T,V);
V = V(1:9e4);
V0 = V0(1:9e4);
dV = [NaN; diff(V)];
dV0 = [NaN; diff(V0)];


% set up the neuron

x = makeSTGNeuron;

% set up the fitter
p = xfit;
p.x = x;

p.parameter_names = [x.find('*gbar');'Axon.len'; 'CellBody.len'; 'CellBody.radius'; 'CellBody.CalciumMech.f'; 'CellBody.CalciumMech.tau_Ca'];
%      Axon                Soma  
%      A   Kd    L    NaV   A   CaS  CaT  H    KCa  Kd   L
p.lb = [10,  100, 0,  100, 10,  1,   1,   1,   10,  10,  0,  .1, .01, .01, 1,  10];
p.ub = [1e3, 1e3, 10, 3e3, 1e3, 100, 100,  10, 1e3, 2e3, 10,  5,  1,  1,  20, 500];

metrics =  xtools.V2metrics(V0,'spike_threshold',-30,'sampling_rate',10);

% prepare data to be sent to cost function
data = struct;
data.max_V = max(V0);
data.min_V = min(V0);
data.burst_period = metrics.burst_period;
data.duty_cycle = metrics.duty_cycle_mean;
data.V0 = V0;
p.data = data;

% debug
p.options.UseParallel = true;

p.sim_func = @metricsCost;



% run forever, churning out neurons
N = 1e4;
all_params = NaN(N,length(p.ub));
all_cost = NaN(N,1);

file_name = 'PD_models.mat';

if exist(file_name)
	load(file_name)
	start_idx = find(isnan(all_cost),1,'first');
else
	start_idx = 1;
end

p.options.MaxTime = MaxTime;
p.options.Display = 'iter';


for i = start_idx:N
	disp(['Starting with random seed #' strlib.oval(i)])
	try
		p.seed = rand(length(p.seed),1).*(p.ub - p.lb) + p.lb;
		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(p.parameter_names,p.seed);
		this_cost = p.sim_func(x, p.data);

		if this_cost < .5


			disp(['Final Cost for this seed is ' strlib.oval(all_cost(i))])

			all_params(i,:) = p.seed;
			all_cost(i) = this_cost;

			save(file_name,'all_params','all_cost')
		else
			disp('Cost too high, skipping...')
		end

	catch
		disp('Something went wrong. Ouch. ')
	end


end
