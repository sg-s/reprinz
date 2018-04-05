% test script for matlab wrapper 

% sets up a single Prinz-like neuron

% conversion from Prinz to phi
A = 0.0628;
vol = A; 
tau_Ca = 200;
phi = .1;

x = xolotl;
x.cleanup;
x.add('AB','compartment','V',-65,'Ca',0.02,'Cm',10,'A',A,'vol',vol,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tau_Ca',tau_Ca);
x.AB.add('prinz/NaV','gbar',1000,'E',50);
x.AB.add('prinz/CaT','gbar',25,'E',30);
x.AB.add('prinz/CaS','gbar',60,'E',30);
x.AB.add('prinz/ACurrent','gbar',500,'E',-80);
x.AB.add('prinz/KCa','gbar',50,'E',-80);
x.AB.add('prinz/Kd','gbar',1000,'E',-80);
x.AB.add('prinz/HCurrent','gbar',.1,'E',-20);

x.t_end = 25e3;
x.dt = .1;
x.sim_dt = .1;

x.transpile;
x.compile;
