
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
x.skip_hash = true;

	% x.add('Soma','compartment','radius',r_soma,'len',L_soma,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tree_idx',0);
	% x.add('Neurite','compartment','radius',r_neurite,'len',L_neurite,'phi',phi,'Ca_out',3000,'Ca_in',0.05);

	% make another copy that we will integrate using the old method
	x.add('Soma','compartment','radius',r_soma,'len',L_soma,'phi',phi,'Ca_out',3000,'Ca_in',0.05);
	x.add('Neurite','compartment','radius',r_neurite,'len',L_neurite,'phi',phi,'Ca_out',3000,'Ca_in',0.05);

	prefix = 'prinz-approx/';
	channels = {'ACurrent','CaS','CaT','HCurrent','KCa','Kd','NaV'};
	g =           [94;      11;    2 ;    .1;      10;  200;  600];
	E =           [-80;      30;    30;  -20;     -80;   -80;  50 ];

	compartments = x.find('compartment');
	for i = 1:length(channels)
		for j = 1:length(compartments)
			x.(compartments{j}).add([prefix channels{i}],'gbar',g(i),'E',E(i));
		end
	end

	x.Soma.NaV.gbar = 0;

	x.connect('Neurite','Soma',axial_resitivity);

x.skip_hash = false; x.sha1hash;


x.dt = .1;
x.sim_dt = .1;
x.t_end = 9e3;

x.synapses(2).gbar = 500;
x.synapses(1).gbar = 100;
