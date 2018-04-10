
x = make2C;
x.transpile;
x.compile;
x.skip_hash_check = true;
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = x.find('Neurite*gbar');

M = length(p.parameter_names);

seed = x.get(p.parameter_names);


% neuron conductances
ub = [500; 100; 100; .5; 100; 1250; 2000];
lb = ub*0 + 1e-2;


p.seed = rand(M,1).*ub; % random seed
x.set(p.parameter_names,p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @two_comp_cost_func;


N = 1e3;
n_epochs = 1;
all_g = NaN(M,N);
all_cost = NaN(N,1);
cost_wo_NaV = NaN(N,1);
spike_height = NaN(N,1);

file_name = ['reprinz_2c_with_spike_ht_' getComputerName '.mat'];

if exist(file_name)
	load(file_name)
	start_idx = find(isnan(all_cost),1,'first');
else
	start_idx = 1;
end

p.options.MaxTime = 300;
p.options.Display = 'none';

for i = start_idx:N
	disp(['Starting with random seed #' oval(i)])
	try
		p.seed = rand(M,1).*ub; % random seed
		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(p.parameter_names,p.seed);	
		all_cost(i) = p.sim_func(x);


		disp(['Final Cost for this seed is ' oval(all_cost(i))])

		all_g(:,i) = p.seed;

		x.set(x.find('Neurite*gbar'),all_g(:,i));	
		x.set(x.find('Soma*gbar'),all_g(:,i));
		x.Soma.NaV.gbar = 0;
		[cost_wo_NaV(i),~,spike_height(i)] = two_comp_cost_func(x);
		 

		save(file_name,'all_g','all_cost','cost_wo_NaV','spike_height')

	catch
		disp('Something went wrong. Ouch. ')
	end



end
