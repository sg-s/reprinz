% test script for matlab wrapper

% sets up a single Prinz-like neuron

function x =  make_neuron_liu()

% conversion from Prinz to phi
% assume spherical geometry
r = 0.0707;
vol = (4/3)*pi*(r^3);
A = 4*pi*r*r;

x = xolotl;
x.add('compartment', 'AB','V',-65,'Ca',0.02,'Cm',10,'A',A,'vol',vol,'Ca_out',3000);
x.AB.add('liu-approx/NaV','gbar',1000,'E',30);
x.AB.add('liu-approx/CaT','gbar',25,'E',30);
x.AB.add('liu-approx/CaS','gbar',60,'E',30);
x.AB.add('liu-approx/ACurrent','gbar',500,'E',-80);
x.AB.add('liu-approx/KCa','gbar',50,'E',-80);
x.AB.add('liu-approx/Kd','gbar',1000,'E',-80);
x.AB.add('liu-approx/HCurrent','gbar',.1,'E',-20);
x.AB.add('Leak','gbar',.1,'E',-50);
x.AB.add('CalciumMech1');

x.t_end = 20e3;
x.dt = .1;
x.sim_dt = .1;
