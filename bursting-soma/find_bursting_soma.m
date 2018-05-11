
x = make_bursting_soma();
x.sha1hash;
p = procrustes('particleswarm');
p.x = x;

p.parameter_names = {'Soma.ACurrent.gbar', 'Soma.CaS.gbar', 'Soma.CaT.gbar', 'Soma.HCurrent.gbar', 'Soma.KCa.gbar' , 'Soma.Kd.gbar' , 'Soma.Leak.gbar' };

M = length(p.parameter_names);

seed = x.get(p.parameter_names);


% neuron conductances
%      A   CaS  CaT   H    KCa  Kd    Leak    
ub = [1e3  600  200  10   2e3   1e3   10    ];
lb = [0    0    0    0    0     0     0     ];


p.seed = rand(M,1).*ub(:); % random seed
x.set(p.parameter_names,p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @bursting_soma_cost_func;


N = 1e4;
n_epochs = 3;
all_g = NaN(M,N);
all_cost = NaN(N,1);

file_name = ['bursting_soma_' getComputerName '.mat'];

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
		p.seed = rand(M,1).*ub(:); % random seed
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
