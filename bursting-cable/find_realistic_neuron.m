
clear x
x = make_realistic_cable();

p = procrustes('particleswarm');
p.x = x;

p.parameter_names = {'CellBody.ACurrent.gbar', 'CellBody.CaS.gbar', 'CellBody.CaT.gbar', 'CellBody.HCurrent.gbar', 'CellBody.KCa.gbar' , 'CellBody.Kd.gbar' , 'CellBody.NaV.gbar' };

M = length(p.parameter_names);

seed = x.get(p.parameter_names);


% neuron conductances
%    A     CaS  CaT   H   KCa  Kd    NaV  
ub = [500  100   50    1  1e3  1e3   2e3 ];
lb = [1    10    1    .1  1    1     1   ];


p.seed = rand(M,1).*ub(:); % random seed
x.set(p.parameter_names,p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @realistic_cost_func;


N = 1e4;
n_epochs = 3;
all_g = NaN(M,N);
all_cost = NaN(N,1);

file_name = ['realistic_neuron_' getComputerName '.mat'];

if exist(file_name)
	load(file_name)
	start_idx = find(isnan(all_cost),1,'first');
else
	start_idx = 1;
end

p.options.MaxTime = 600;
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
