


% first, we create our xolotl object
x = xolotl.examples.networks.pyloric;
p = xfit;
p.x = x;

p.options.UseParallel = true;



% we assign a cost function
p.SimFcn = @STG_cost_function;



% we optimize over all maximal conductances and synapses
p.FitParameters = [x.find('*gbar'); x.find('*gmax')];

% bounds
seed = x.get(p.FitParameters);
ub = 0*seed;
lb = 0*seed;

% neuron conductances
%                  A    CaS   CaT  H  KCa  Kd    L   NaV
ub(1:24) = repmat([500; 100; 100; .5; 100; 1250; 1 ; 4000],3,1);
lb(1:24) = repmat([100; 0  ; 0  ; 0 ;   0;  250; 0 ; 400 ],3,1);


% synapses
ub(25:31) = 100; % nS
lb(25:31) = 0; % nS


p.lb = lb;
p.ub = ub;

p.SaveParameters = p.FitParameters;

p.SaveWhenCostBelow = 1;






