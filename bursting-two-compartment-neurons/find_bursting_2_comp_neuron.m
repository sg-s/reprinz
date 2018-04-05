
clear all
x = make2C;
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = x.find('Neurite*gbar');

seed = x.get(p.parameter_names);


% neuron conductances
ub = [500; 100; 100; .5; 100; 1250; 2000];
lb = ub*0 + 1e-2;


p.seed = rand(length(p.parameter_names),1).*ub; % random seed
x.set(p.parameter_names,p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @two_comp_cost_func;


N = 1e3;
n_epochs = 3;
all_g = NaN(28,N);
all_metrics = NaN(18,N);
all_cost = NaN(N,1);

file_name = ['reprinz_2c_' getComputerName '.mat'];

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
