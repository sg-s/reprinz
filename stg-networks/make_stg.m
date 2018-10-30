% test script for matlab wrapper 

% this sets up the STG network 
% as in Fig 2e of this paper:
% Prinz ... Marder Nat Neuro 2004
% http://www.nature.com/neuro/journal/v7/n12/abs/nn1352.html

function x = make_STG(prefix)

% conversion from Prinz to phi
A = 0.0628;

channels = {'NaV','CaT','CaS','ACurrent','KCa','Kd','HCurrent','Leak'};

gbar(:,1) = [1000 25  60 500  50  1000 .1  0 ];
gbar(:,2) = [1000 0   40 200  0   250  .5 .3];
gbar(:,3) = [1000 24  20 500  0   1250 .5 .1];
E =         [50   30  30 -80 -80 -80   -20 -50];

x = xolotl;

x.add('compartment','AB','Cm',10,'A',A);
x.add('compartment','LP','Cm',10,'A',A);
x.add('compartment','PY','Cm',10,'A',A);

compartments = x.find('compartment');
for j = 1:length(compartments)

	% add Calcium mechanism
	x.(compartments{j}).add('CalciumMech1');

	for i = 1:length(channels)-1

		x.(compartments{j}).add([prefix channels{i}],'gbar',gbar(i,j),'E',E(i));
	end
	x.(compartments{j}).add(channels{8},'gbar',gbar(8,j),'E',E(8));
end



% set up synapses as in Fig. 2e
x.connect('AB','LP','Chol','gbar',30);
x.connect('AB','PY','Chol','gbar',3);
x.connect('AB','LP','Glut','gbar',30);
x.connect('AB','PY','Glut','gbar',10);
x.connect('LP','PY','Glut','gbar',1);
x.connect('PY','LP','Glut','gbar',30);
x.connect('LP','AB','Glut','gbar',30);


x.t_end = 5e3;


x.t_end = 25e3;
x.dt = .1;
x.sim_dt = .1;
