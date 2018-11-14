% the point of this script
% is to find conditional bursters
% using the neuroDB project
%
% we want to find cells that are spiking tonically
% w/o any synaptic input, but fire in a HCO
% when forced by a conditional bursters


if ~exist('synapse_type','var')
	error('No synapse_type defiend')
end

if exist('bursters.mat')
	load('bursters.mat')
else

	if ~exist('n','var')
		disp('Connecting to neuron DB...')
		n = neuroDB;
		n.prefix = 'prinz/';
	end


	show_these = find(n.results.burst_period > 800 ...
					& n.results.burst_period < 2e3 ...
					& n.results.duty_cycle_mean > .3 ...
					& n.results.duty_cycle_mean < .5 ...
					& n.results.burst_period_std < .1 ...
		            & n.results.duty_cycle_std./n.results.duty_cycle_mean < .1 ...
		            & n.results.min_V_in_burst_mean > n.results.min_V_mean ...
		            & n.results.min_V_mean < -60 ...
		            & n.results.spike_peak_std./n.results.spike_peak_mean < .1); 

	all_g = n.results(show_these,:);

end




x = make_cb_network(synapse_type);
mkdir(synapse_type)
p = procrustes('particleswarm');

p.x = x;

p.parameter_names = x.find('CB*gbar');


p.options.UseParallel = true;
seed = x.get(x.find('*gbar'));
ub = 0*seed;
lb = 0*seed;


%      A    CaS   CaT  H  KCa  Kd    L   NaV
ub = [500; 100; 100; .5; 100; 1250; 1 ; 4000];
lb = [100; 0  ; 0  ; 0 ;   0;  250; 0 ; 400 ];



p.seed = rand(8,1).*ub; % random seed
p.lb = lb;
p.ub = ub;

p.sim_func = @conditional_burster_cf;

n_epochs = 2;


p.options.MaxTime = 500;
p.options.Display = 'iter';


while true

	try

		g = all_g(randi(length(all_g)),:);
		p.x.set('Burster*gbar',g);
		p.seed = rand(1,8)*(ub- lb) + lb;

		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(p.parameter_names,p.seed);	

		[C, Ca_burster, Ca_CB0, Ca_CB1] = p.sim_func(p.x);

		if C < 10
			% save it
			disp('Found one!')

			% measure the Calcium levels in the two cells
			% with and without the synapse 


			file_name = [synapse_type '/' GetMD5(p.seed) '.mat'];
			g = x.get('*gbar');
			save(file_name,'g','C','Ca_burster','Ca_CB0','Ca_CB1','-v7.3','-nocompression')

		end

	catch
		disp('Something went wrong. Ouch. ')
	end

end
