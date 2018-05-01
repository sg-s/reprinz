
make_stg;
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


N = 1e3;
n_epochs = 5;
all_g = NaN(28,N);
all_metrics = NaN(18,N);
all_cost = NaN(N,1);

file_name = ['reprinz_' getComputerName '.mat'];

if exist(file_name)
	load(file_name)
	start_idx = find(isnan(all_cost),1,'first');
else
	start_idx = 1;
end

p.options.MaxTime = 300;
p.options.Display = 'iter';

for i = start_idx:N
	disp(['Starting with random seed #' oval(i)])
	try
		p.seed = rand(28,1).*ub;
		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(x.find('*gbar'),p.seed);	
		[all_cost(i),~,all_metrics(:,i)] = p.sim_func(x);


		disp(['Final Cost for this seed is ' oval(all_cost(i))])

		all_g(:,i) = p.seed;
		 

		save(file_name,'all_g','all_metrics','all_cost')

	catch
		disp('Something went wrong. Ouch. ')
	end



end
