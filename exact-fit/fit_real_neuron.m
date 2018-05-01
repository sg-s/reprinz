
load('isolated_PD')
V = V(1:9e4);


g_axial = 10; % nS

% geometry of soma is a sphere
r = 0.0707;
A = 4*pi*r*r;
vol = (4/3)*pi*r*r*r;
f = 14.96; % uM/nA
tau_Ca = 200;
F = 96485; % Faraday constant in SI units
phi = (2*f*F*vol)/tau_Ca;

x = xolotl;
x.add('Soma','compartment','Cm',10,'A',A,'vol',vol,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tau_Ca',tau_Ca);

x.add('Neurite','compartment','Cm',10,'A',A/10,'vol',vol/10,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tau_Ca',tau_Ca);

x.connect('Neurite','Soma',g_axial);

E_K = -80;

g0 = 1e-1*rand(7,1);
x.Soma.add('prinz-approx/NaV','gbar',g0(1),'E',50);
x.Soma.add('prinz-approx/CaT','gbar',g0(2),'E',30);
x.Soma.add('prinz-approx/CaS','gbar',g0(3),'E',30);
x.Soma.add('prinz-approx/ACurrent','gbar',g0(4),'E',E_K);
x.Soma.add('prinz-approx/KCa','gbar',g0(5),'E',E_K);
x.Soma.add('prinz-approx/Kd','gbar',g0(6),'E',E_K);
x.Soma.add('prinz-approx/HCurrent','gbar',g0(7),'E',-20);
x.Soma.add('Leak','gbar',1,'E',-50);

% now add some conductances and controllers to the neurite
g0 = 1e-1*rand(7,1);
x.Neurite.add('prinz-approx/NaV','gbar',g0(1),'E',50);
x.Neurite.add('prinz-approx/CaT','gbar',g0(2),'E',30);
x.Neurite.add('prinz-approx/CaS','gbar',g0(3),'E',30);
x.Neurite.add('prinz-approx/ACurrent','gbar',g0(4),'E',E_K);
x.Neurite.add('prinz-approx/KCa','gbar',g0(5),'E',E_K);
x.Neurite.add('prinz-approx/Kd','gbar',g0(6),'E',E_K);
x.Neurite.add('prinz-approx/HCurrent','gbar',g0(7),'E',-20);
x.Neurite.add('Leak','gbar',.099,'E',-50);

x.t_end = 15e3;


x.dt = .1;
x.sha1hash
p = procrustes('particleswarm');
p.x = x;
p.V_clamp = [V];

p.parameter_names = [x.find('Soma*gbar'); 'Soma.tau_Ca'];

M = length(p.parameter_names);


% neuron conductances
p.ub = [1e3  1e3  1e3  .5  100  3e4   100 3e3  300];
p.lb = [.1   .1   .1   .1   .1   500  10  100  20];


p.seed = rand(M,1).*p.ub(:); % random seed

p.sim_func = @real_data_cost_func;
p.options.MaxTime = 300;




