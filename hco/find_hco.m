% the point of this script
% is to find half-center oscilaltor networks
% using the neuroDB project
%



if ~exist('synapse_type','var')
	error('No synapse_type defiend')
end

cb_db = dir(['../ConditionalBurster/' synapse_type '/*.data']);
assert(length(cb_db) == 1, 'Expected one .data file in ../ConditionalBurster/..')

cb_db = Data([cb_db.folder filesep cb_db.name]);
cb_db = cb_db.filter(cb_db.C == 0);



x = makeHCONetwork(synapse_type);
mkdir(synapse_type)
p = procrustes('particleswarm');

p.x = x;

p.parameter_names = x.find('Cell1*gbar');


p.options.UseParallel = true;



seed = x.get(x.find('Cell1*gbar'));
ub = 0*seed;
lb = 0*seed;


%      A    CaS   CaT  H  KCa  Kd    L   NaV
ub = [500; 100; 100; .5; 100; 1250; 1 ; 4000];
lb = [100; 0  ; 0  ; 0 ;   0;  250; 0 ; 400 ];



p.seed = rand(8,1).*ub; % random seed
p.lb = lb;
p.ub = ub;

p.sim_func = @hco_cost_func;

n_epochs = 1;


p.options.MaxTime = 500;
p.options.Display = 'iter';


while true

	try

		g_cb = cb_db.g(randi(cb_db.size),:);
		g_cb(1:8) = [];
		g_cb(end) = [];

		p.x.set('Cell1*gbar',g_cb);
		p.seed = rand(1,8)*(ub- lb) + lb;

		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(p.parameter_names,p.seed);	

		[C, metrics] = p.sim_func(p.x);

		if C == 0 
			% save it
			disp('Found one!')

			% measure the Calcium levels in the two cells
			% with and without the synapse 

			H = GetMD5(which(mfilename),'File');
			g = x.get('*gbar');
			file_name = [synapse_type filesep H '_' GetMD5(g) '.mat'];
			
			metrics.g = g;

			metrics = Data(metrics);

			metrics.save;

		end

	catch
		disp('Something went wrong. Ouch. ')
	end

end
