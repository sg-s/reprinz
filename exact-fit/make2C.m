
function x = make2C()


axial_resitivity = .001; % MOhm mm; 

% geometry is a cylinder
r_neurite = .01;
L_neurite = .35; % mm, from Otopalik et al
r_soma = .025;
L_soma = .05; % mm, Jason's guesstimates
f = 14.96; % uM/nA
phi = 100;

x = xolotl;
% x.skip_hash = true;


	% make compartments 
	x.add('compartment','CellBody','radius',r_soma,'len',L_soma,'phi',phi,'Ca_out',3000,'Ca_in',0.05);
	x.add('compartment','Neurite','radius',r_neurite,'len',L_neurite,'phi',phi,'Ca_out',3000,'Ca_in',0.05);


	x.CellBody.add('prinz-temperature/ACurrent','gbar',10);
	x.CellBody.add('prinz-temperature/CaS','gbar',10);
	x.CellBody.add('prinz-temperature/CaT','gbar',10);
	x.CellBody.add('liu-temperature/HCurrent','gbar',10);
	x.CellBody.add('prinz-temperature/KCa','gbar',10);
	x.CellBody.add('prinz-temperature/Kd','gbar',10);
	x.CellBody.add('Leak','gbar',10);


	x.Neurite.add('prinz-temperature/ACurrent','gbar',10);
	x.Neurite.add('prinz-temperature/CaS','gbar',10);
	x.Neurite.add('prinz-temperature/CaT','gbar',10);
	x.Neurite.add('liu-temperature/HCurrent','gbar',10);
	x.Neurite.add('prinz-temperature/KCa','gbar',10);
	x.Neurite.add('prinz-temperature/Kd','gbar',10);
	x.Neurite.add('prinz-temperature/NaV','gbar',10);
	x.Neurite.add('Leak','gbar',10);

	x.CellBody.tree_idx = 0;
	x.Neurite.tree_idx = 1; 
	x.connect('Neurite','CellBody',axial_resitivity);

	x.temperature_ref = 22;


% x.skip_hash = false; x.md5hash;


x.dt = .1;
x.sim_dt = .1;
x.t_end = 9e3;
