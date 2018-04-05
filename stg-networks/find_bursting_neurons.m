
% lazy code, we're going to reuse the STG network code
% and turn off all the synapses

make_neuron;
p = procrustes('particleswarm');

x.AB.phi = 1;

p.x = x;


p.parameter_names = [x.find('*gbar'); 'AB.Ca_in' ; 'AB.Ca_out'] ;

seed = x.get(p.parameter_names);


% neuron conductances
ub = [500; 100; 100; .5; 100; 2e3; 2e3; .1; 5000];
lb(1:7) = 1e-2;
lb(8) = .01; % Ca_in
lb(9) = .2; % Ca_out

p.seed = rand(9,1).*ub; % random seed
x.set(p.parameter_names,p.seed);
p.lb = lb;
p.ub = ub;

p.sim_func = @bursting_neuron;