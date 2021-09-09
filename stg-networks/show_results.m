x = xolotl.examples.networks.pyloric;

load prinz_networks.mat	
alldata = structlib.purge(alldata,alldata.cost~=0);


for i = 1:length(alldata.cost)

	

	x.set(x.find('*gbar'),alldata.gbar(i,:));	
	x.set(x.find('*gmax'),alldata.gmax(i,:));	
	x.plot

	drawnow
	pause(3)
	close all
	

end
