% this cost function can be used to fit a two compartment 
% model to data using metrics gleaned from the data


function C = metricsCost(x, data)

x.reset;

x.closed_loop = true;
x.dt = .1;
x.t_end = 10e3;
x.integrate;

V = x.integrate;

if any(isnan(V(:)))
	C = 100;
	return
end

C = 0;

% measure metrics in Axon
metrics = xtools.V2metrics(V(:,1),'sampling_rate',10);

C = C + xfit.binCost([data.burst_period - 100, data.burst_period + 100],metrics.burst_period);
C = C + xfit.binCost([data.duty_cycle - .05, data.duty_cycle + .05],metrics.duty_cycle_mean);

C = C + 5*xfit.binCost([data.n_spikes_per_burst - 1, data.n_spikes_per_burst + 3],metrics.n_spikes_per_burst_mean);


% measure minimum and maximum in soma
C = C + 5*xfit.binCost([data.min_V - 2, data.min_V + 2],min(V(:,2)));
C = C + 5*xfit.binCost([data.max_V - 2, data.max_V + 2],max(V(:,2)));

% also compare voltage traces directly
C = C + 10*xtools.voltageCost(data.V0,V(:,2),100);

% also measure the minimum voltage b/w spikes on the PD
spiketimes = xtools.findNSpikeTimes(V(:,1),xtools.findNSpikes(V(:,1)));

spiketimes(spiketimes > length(V)) = [];

% burst cost
BC = 0;


if length(spiketimes) > 2

	for i = 2:length(spiketimes)


		if (spiketimes(i) - spiketimes(i-1))*x.dt*1e-3 < .3
			% this is in a burst
			min_V = min(V(spiketimes(i-1):spiketimes(i),2));
			BC = BC + xfit.binCost(data.V_bw_spikes_range, min_V);

		end
	end
	BC = BC/length(spiketimes);
end



% spike peaks
BC = BC + xfit.binCost(data.spike_peaks,max(V(:,2)));

C = C + BC;