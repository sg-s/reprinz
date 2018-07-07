function x =  make_neuron()


x = xolotl;
x.add('compartment','AB');
x.AB.add('CalciumMech1');
x.AB.add('prinz/NaV','gbar',1000,'E',50);
x.AB.add('prinz/CaT','gbar',25,'E',30);
x.AB.add('prinz/CaS','gbar',60,'E',30);
x.AB.add('prinz/ACurrent','gbar',500,'E',-80);
x.AB.add('prinz/KCa','gbar',50,'E',-80);
x.AB.add('prinz/Kd','gbar',1000,'E',-80);
x.AB.add('liu/HCurrent','gbar',.1,'E',-20);
x.AB.add('Leak','E',-50);

x.t_end = 20e3;
x.dt = .1;
x.sim_dt = .1;
