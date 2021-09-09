% makes a realistic cable with some options
function x = make_realistic_cable(varargin)

% options and defaults
options.N = 10;
options.n_channels = 7;
options.uniform_nav = false;
options.start_axon = 7;

options = corelib.parseNameValueArguments(options,varargin{:});


axial_resitivity = 1e-3; % MOhm mm; 

% geometry is a cylinder
r_neurite = .01; % default value
L_neurite = .35*5; % mm, from Otopalik et al
r_soma = .025;
L_soma = .05; % mm, Jason's guesstimates
phi = 1;

shell_thickness = .01; % 10 micrometres

x = xolotl;
% make 3 neurons
x.add('compartment','CellBody','radius',r_soma,'len',L_soma,'Ca_out',3000,'tree_idx',0,'shell_thickness',shell_thickness);

x.CellBody.add('buchholtz/CalciumMech','phi',phi);

if options.N > 1
	x.add('compartment','Neurite','radius',r_neurite,'len',L_neurite,'Ca_out',3000,'shell_thickness',shell_thickness);
	x.Neurite.add('buchholtz/CalciumMech','phi',phi);
end


prefix = 'prinz/';
channels = {'ACurrent','CaS',  'CaT','HCurrent','KCa','Kd','NaV'};
%g =           [104;     11.76;  4.7 ;   .1;      490;  999;  1726];
g =           [104;     11.76;  4.7 ;   .1;      390;  250;  2e3];
E =           [-80;      30;   30;      -20;     -80;   -80;  50 ];

compartments = x.find('compartment');
for i = 1:options.n_channels
	for j = 1:length(compartments)
		x.(compartments{j}).add([prefix channels{i}],'gbar',g(i),'E',E(i));

	end
end

if options.N > 1
	x.slice('Neurite',options.N);

	comp_names = x.find('compartment');
	if ~options.uniform_nav
		for i = 1:options.start_axon
			try
				x.(comp_names{i}).NaV.gbar = 0;
			catch
			end
		end
	end

	x.connect('Neurite','CellBody');
end



x.dt = .1;
x.sim_dt = .1;