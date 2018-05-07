% test script for matlab wrapper

% sets up a single Prinz-like neuron

function x =  make_neuron_liu()

% conversion from Prinz to phi
% assume spherical geometry
r = 0.0707;
vol = (4/3)*pi*(r^3);
A = 4*pi*r*r;

f = 14.96; % uM/nA
tau_Ca = 200;
F = 96485; % Faraday constant in SI units
phi = (2*f*F*vol)/tau_Ca;

x = xolotl;
x.add('AB','compartment','V',-65,'Ca',0.02,'Cm',10,'A',A,'vol',vol,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tau_Ca',tau_Ca);
x.AB.add('liu-approx/NaV','gbar',1000,'E',30);
x.AB.add('liu-approx/CaT','gbar',25,'E',30);
x.AB.add('liu-approx/CaS','gbar',60,'E',30);
x.AB.add('liu-approx/ACurrent','gbar',500,'E',-80);
x.AB.add('liu-approx/KCa','gbar',50,'E',-80);
x.AB.add('liu-approx/Kd','gbar',1000,'E',-80);
x.AB.add('liu-approx/HCurrent','gbar',.1,'E',-20);
x.AB.add('Leak','gbar',.1,'E',-50);

x.t_end = 20e3;
x.dt = .1;
x.sim_dt = .1;
