
load('isolated_PD')
V_clamp = V(1:9e4);
%V_clamp = [NaN*V_clamp V_clamp ];
x = make2C;

x.t_end = 9e3;
x.dt = .1;
x.sim_dt = .1;
x.V_clamp = repmat(V_clamp,1,2);
x.V_clamp(:,2) = NaN;

p = procrustes('particleswarm');
p.x = x;

p.parameter_names = [x.find('Soma*gbar'); 'Soma.tau_Ca'; 'Soma.Ca'];

M = length(p.parameter_names);


% neuron conductances
%       A     CaS  CaT  H   KCa  Kd    NaV   tau_Ca  Ca 
p.ub = [1e3  1e3  1e3  .5  100   2e3   2e3   500     3e3  ];
p.lb = [.1   .1   .1   .1   .1   100   500   20      .05  ];


p.seed = rand(M,1).*p.ub(:); % random seed

p.sim_func = @exact_fit_cost_func;
p.options.MaxTime = 300;




