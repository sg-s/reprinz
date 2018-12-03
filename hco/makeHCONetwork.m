% makes a half-center oscillator network
% with two neurons and two synapses connecting them
% the two neurons and the two synapses are identical 
% 
function x = makeHCONetwork(synapse_type)


% conversion from Prinz to phi
A = 0.0628;

channels = {'NaV','CaT','CaS','ACurrent','KCa','Kd','HCurrent','Leak'};

gbar(:,1) = [1000 25  60 500  50  1000 .1  0 ];
gbar(:,2) = [1000 0   40 200  0   250  .5 .3];

E =         [50   30  30 -80 -80 -80   -20 -50];

x = xolotl;

x.add('compartment','Cell1','Cm',10,'A',A);
x.add('compartment','Cell2','Cm',10,'A',A);

prefix = 'prinz/';

compartments = x.find('compartment');
for j = 1:length(compartments)

	% add Calcium mechanism
	x.(compartments{j}).add('prinz/CalciumMech');

	for i = 1:length(channels)-1

		x.(compartments{j}).add([prefix channels{i}],'gbar',gbar(i,j),'E',E(i));
	end
	x.(compartments{j}).add(channels{8},'gbar',gbar(8,j),'E',E(8));
end


% connect cells to each other 
x.connect('Cell1','Cell2',synapse_type);
x.connect('Cell2','Cell1',synapse_type);