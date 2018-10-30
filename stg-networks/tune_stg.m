

if ~exist('prefix','var')
	error('No prefix defiend')
end

x = make_stg(prefix);
mkdir(prefix)
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = x.find('*gbar');

p.options.UseParallel = true;

seed = x.get(x.find('*gbar'));
ub = 0*seed;
lb = 0*seed;

% neuron conductances
%                  A    CaS   CaT  H  KCa  Kd    L   NaV
ub(1:24) = repmat([500; 100; 100; .5; 100; 1250; 1 ; 4000],3,1);
lb(1:24) = repmat([100; 0  ; 0  ; 0 ;   0;  250; 0 ; 400 ],3,1);

% close to reference network
%                  A    CaS   CaT  H   KCa  Kd    L   NaV
ub(1:24) = repmat([500; 60;  25;   .5;  50; 1250; .1; 1000],3,1);
lb(1:24) = repmat([200; 20 ; 0 ;  .1 ;   0;  250; 0 ; 1e3 ],3,1);

% synapses
ub(25:31) = 100; % nS
lb(25:31) = 0; % nS



p.seed = rand(31,1).*ub; % random seed
x.set(x.find('*gbar'),p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @STG_cost_function;

n_epochs = 5;

if exist(file_name)
	load(file_name)
	start_idx = find(isnan(all_cost),1,'first');
else
	start_idx = 1;
end

p.options.MaxTime = 300;
p.options.Display = 'iter';

idx = 1;

while true
	disp(['Starting with random seed #' oval(idx)])
	try
		p.seed = rand(31,1).*(ub-lb) + lb;
		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(x.find('*gbar'),p.seed);	
		[all_cost,~,all_metrics] = p.sim_func(x);


		disp(['Final Cost for this seed is ' oval(all_cost)])

		all_g = p.seed;

		file_name = [prefix GetMD5(p.seed) '.mat'];

		save(file_name,'all_g','all_metrics','all_cost','-v7.3','-nocompression')

		idx = idx + 1;

	catch
		disp('Something went wrong. Ouch. ')
	end

end
