
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

p.seed = rand(length(p.seed),1).*(p.ub - p.lb) + p.lb;

