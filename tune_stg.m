
make_stg
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = x.find('*gbar');

seed = x.get(x.find('*gbar'));

p.seed = []; % random seed
p.lb = seed/10;
p.ub = seed*5;

p.sim_func = @STG_cost_function;

