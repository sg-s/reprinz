% makes a realistic cable with some options
function x = make2C(varargin)



axial_resitivity = 1e-3; % MOhm mm; 

r_soma = .025;
L_soma = .05; % mm, Jason's guesstimates
phi = 1;

r_neurite = .01;
L_neurite = .35; % mm, from Otopalik et al

shell_thickness = .01; % 10 micrometres

x = xolotl;

x.add('compartment','CellBody','radius',r_soma,'len',L_soma,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tree_idx',0,'shell_thickness',shell_thickness);

x.add('compartment','Neurite','radius',r_neurite,'len',L_neurite,'phi',phi,'Ca_out',3000,'Ca_in',0.05);



% add the channels
x.CellBody.add('prinz/ACurrent','E',-80,'gbar',1);
x.CellBody.add('prinz/CaS','E',0,'gbar',1);
x.CellBody.add('prinz/CaT','E',0,'gbar',1);
x.CellBody.add('liu/HCurrent','E',-20,'gbar',1);
x.CellBody.add('prinz/KCa','E',-80,'gbar',1);
x.CellBody.add('prinz/Kd','E',-80,'gbar',1);
x.CellBody.add('Leak','E',-50,'gbar',1);

% add the channels
x.Neurite.add('prinz/ACurrent','E',-80,'gbar',1);
x.Neurite.add('prinz/CaS','E',0,'gbar',1);
x.Neurite.add('prinz/CaT','E',0,'gbar',1);
x.Neurite.add('liu/HCurrent','E',-20,'gbar',1);
x.Neurite.add('prinz/KCa','E',-80,'gbar',1);
x.Neurite.add('prinz/Kd','E',-80,'gbar',1);
x.Neurite.add('Leak','E',-50,'gbar',1);
x.Neurite.add('prinz/NaV','E',50,'gbar',1e3);

x.CellBody.tree_idx = 0;
x.connect('CellBody','Neurite');


x.dt = 1;
x.sim_dt = 1;
x.t_end = 15e3;

x.md5hash;