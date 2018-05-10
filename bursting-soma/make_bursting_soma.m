% this function makes a soma compartment 
% that will eventually be part of a mulit-compartment model
% with NaV, etc. 

function x = make_bursting_soma()


r_soma = .025;
L_soma = .05; % mm, Jason's guesstimates
phi = 1;

shell_thickness = .01; % 10 micrometres

x = xolotl;

x.add('Soma','compartment','radius',r_soma,'len',L_soma,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tree_idx',0,'shell_thickness',shell_thickness);

prefix = 'liu/';
channels = {'ACurrent','CaS','CaT','HCurrent','KCa','Kd','Leak'};
g =           [500;     100;   25 ;   1;      10;   1e3;  1500];
E =           [-80;      30;   30;   -20;     -80;  -80;  -50 ];

for i = 1:length(channels)-1
	x.Soma.add([prefix channels{i}],'gbar',g(i),'E',E(i));
end
i = i + 1;
x.Soma.add([channels{i}],'gbar',g(i),'E',E(i));


x.dt = .2;
x.sim_dt = .2;
x.t_end = 15e3;

x.sha1hash;