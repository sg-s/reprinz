
x = make_neuron;
x.sha1hash;
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = x.find('*gbar');

M = length(p.parameter_names);

seed = x.get(p.parameter_names);


% neuron conductances
experimental_range = [500; 100; 100; .5; 100; 1250; 2000];
ub = [1e3 1e3 200 10 1e3 2e3 2e3];
lb = ub*0 + 1e-2;


p.seed = rand(M,1).*ub; % random seed
x.set(p.parameter_names,p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @bursting_cost_func;


N = 1e4;
n_epochs = 1;
all_g = NaN(M,N);
all_cost = NaN(N,1);

file_name = ['reprinz_1c_' getComputerName '.mat'];

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
	%try
		p.seed = rand(M,1).*experimental_range; % random seed
		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(p.parameter_names,p.seed);	
		all_cost(i) = p.sim_func(x);


		disp(['Final Cost for this seed is ' oval(all_cost(i))])

		all_g(:,i) = p.seed;
		 

		save(file_name,'all_g','all_cost')

	% catch
	% 	disp('Something went wrong. Ouch. ')
	% end



end
