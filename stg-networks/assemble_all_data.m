%% 


% first gather all the data
system([' scp -r  enkidu.bio.brandeis.edu:~/code/reprinz/reprinz_enkidu.mat ./'])
system([' scp -r  neuromancer.bio.brandeis.edu:~/code/reprinz/reprinz_meuromancer.mat ./'])

allfiles = dir('reprinz*.mat');

g = NaN(28,2e3);
C = NaN(2e3,1);
M = NaN(18,2e3);

for i = 1:length(allfiles)
	load(allfiles(i).name)

	keep_these = ~isnan(all_cost);
	all_cost = all_cost(keep_these);
	all_g = all_g(:,keep_these);
	all_metrics = all_metrics(:,keep_these);

	a = find(isnan(C),1,'first');
	z = a + length(all_cost) - 1;

	% assemble
	C(a:z) = all_cost;
	M(:,a:z) = all_metrics;
	g(:,a:z) = all_g;

end

% trim
a = find(isnan(C),1,'first');
C = C(1:a-1);
M = M(:,1:a-1);
g = g(:,1:a-1);