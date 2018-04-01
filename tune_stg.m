
make_stg
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = x.find('*gbar');

seed = x.get(x.find('*gbar'));
ub = 0*seed;
lb = 0*seed;

% neuron conductances
ub(1:21) = repmat([500; 100; 100; .5; 100; 1250; 2000],3,1);
lb(1:21) = 1e-2;

% synapses
ub(22:28) = 100; % nS
lb(22:28) = 0; % nS



p.seed = rand(28,1).*ub; % random seed
x.set(x.find('*gbar'),p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @STG_cost_function;

