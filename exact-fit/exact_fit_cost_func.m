function C = exact_fit_cost_func(x,data)

inf_cost = 1e3;
C = inf_cost;


% copy params
% x.CellBody.ACurrent.gbar = x.Neurite.ACurrent.gbar; 
% x.CellBody.CaS.gbar = x.Neurite.CaS.gbar;
% x.CellBody.CaT.gbar = x.Neurite.CaT.gbar;
% x.CellBody.HCurrent.gbar = x.Neurite.HCurrent.gbar;
% x.CellBody.KCa.gbar = x.Neurite.KCa.gbar;
% x.CellBody.Kd.gbar = x.Neurite.Kd.gbar;

x.CellBody.Ca = x.Neurite.Ca;
x.CellBody.tau_Ca = x.Neurite.tau_Ca;
x.synapses(2).resistivity = x.synapses(1).resistivity;

x.CellBody.V = -54;
x.Neurite.V = -54;

V = x.integrate;


if nargout == 1

	M = procrustes.V2matrix(V(45000:end,1),[-55 -20],[-0.3 .5]);

	C = 1e4*procrustes.matrixCost(data.LeMassonMatrix,M);


	% also add a cost for the spiketimes 
	N = length(data.spiketimes);
	spiketimes = xolotl.findNSpikeTimes(V(:,1),N,-30);

	for i = 1:length(data.spiketimes)
		if isnan(spiketimes(i))
			this_cost = 20;
		
		else
			this_cost =  abs(spiketimes(i) - data.spiketimes(i));

			if this_cost > 10
				this_cost = 10;
			end

		end

		C = C + this_cost;


	end


	return
end

figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on


subplot(2,1,1); hold on
plot(data.V,'k')
plot(V(:,1),'r')

subplot(2,1,2); hold on
plot(V(:,2))


prettyFig();

