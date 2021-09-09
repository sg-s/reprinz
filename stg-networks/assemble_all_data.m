%% 




allfiles = dir('prinz/*.mat');

N = length(allfiles);

alldata.gbar = NaN(N,24);
alldata.gmax = NaN(N,7);
alldata.metrics = NaN(N,18);
alldata.cost = NaN(N,1);

for i = 1:length(allfiles)

	corelib.textbar(i,length(allfiles))
	load(fullfile(allfiles(i).folder, allfiles(i).name))

	if length(all_cost) > 1
		error('More than one data pt')
	end

	if isnan(all_cost)
		continue
	end

	alldata.gbar(i,:) = all_g(1:24);
	alldata.gmax(i,:) = all_g(25:end);
	alldata.cost(i) = all_cost;
	alldata.metrics(i,:) = all_metrics;

end
