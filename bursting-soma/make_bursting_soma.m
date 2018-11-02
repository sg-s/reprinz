% this function makes a soma compartment
% that will eventually be part of a mulit-compartment model
% with NaV, etc.

function x = make_bursting_soma()


r_soma = .025;
L_soma = .05; % mm, Jason's guesstimates
phi = 1;

shell_thickness = .01; % 10 micrometres

x = xolotl;

x.add('compartment','CellBody','radius',r_soma,'len',L_soma,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tree_idx',0,'shell_thickness',shell_thickness);

% add the channels
x.CellBody.add('prinz/ACurrent','E',-80);
x.CellBody.add('prinz/CaS','E',0);
x.CellBody.add('prinz/CaT','E',0);
x.CellBody.add('liu/HCurrent','E',-20);
x.CellBody.add('prinz/KCa','E',-80);
x.CellBody.add('prinz/Kd','E',-80);
x.CellBody.add('Leak','E',-50);


x.dt = 1;
x.sim_dt = 1;
x.t_end = 15e3;

x.md5hash;
