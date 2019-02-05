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
p = xfit('particleswarm');

p.x = x;

p.parameter_names = x.find('Cell1*gbar');


p.options.UseParallel = true;



seed = x.get(x.find('Cell1*gbar'));
ub = 0*seed;
lb = 0*seed;


%      A   CaS  CaT  H   KCa   Kd   L   NaV
ub = [500; 100; 100; 10; 100; 1250; 1 ; 4000];
lb = [100; 0  ; 0  ; 0 ;   0;  250; 0 ; 400 ];



p.seed = rand(8,1).*(ub-lb) + lb; % random seed
p.lb = lb;
p.ub = ub;

p.sim_func = @hco_cost_func;

n_epochs = 3;


p.options.MaxTime = 500;
p.options.Display = 'iter';


while true

	try

		g_cb = cb_db.g(randi(cb_db.size),:);
		g_cb = g_cb(9:end-1);

		p.x.set('Cell1*gbar',g_cb);
		p.seed = g_cb;

		for j = 1:n_epochs
			p.fit;
		end

		% save
		x.set(p.parameter_names,p.seed);	

		[C, metrics] = p.sim_func(p.x);


		if C < 5
			% save it
			disp('Found one!')

			% measure the Calcium levels in the two cells
			% with and without the synapse 

			H = GetMD5(which(mfilename),'File');
			g = x.get('*gbar');
			file_name = [synapse_type filesep H '_' GetMD5(g) '.mat'];
			
			metrics.g = g;
			metrics.C = C;

			metrics = Data(metrics);

			metrics.save;

		end

	catch
		disp('Something went wrong. Ouch. ')
	end

end
