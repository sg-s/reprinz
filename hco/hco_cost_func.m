% this cost function can be used to find 
% neurons that tonically spike without any inhibitory
% input, and fire bursts between bursts of inhib.
% synaptic input
% 
function [COST, metrics] = hco_cost_func(x,~,~)

COST = 0;

% parameters
A1 = 300; % synapse off firing rate
A2 = 100; % regular spiking cost
B = 10; % burst period
C = 100; % duty cycle
D = 100; % ratio of burt periods for two cells
E = 100; % ratio of firing rates of two cells
F = 200; % phase

x.reset;

% mirror gbars in Cell2
x.set('Cell2*gbar',x.get('Cell1*gbar'));

Ca_burster = NaN;
Ca_CB0 = NaN;
Ca_CB1 = NaN;

% don't change these!

x.dt = .1;
x.sim_dt = .1;

% turn synapse off. now, we want both neurons to be
% tonically spiking. since the two neurons are identical
% it is sifficient to check only one

x.reset;
x.set('*gmax',0)
x.closed_loop = false;

x.t_end = 20e3;
V = x.integrate;
V(1:1e4,:) = [];

m = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);

COST = COST + A1*xfit.binCost([2, 50], m.firing_rate);


COST = COST + A2*xfit.binCost([0, 1e-2],m.isi_std/m.isi_mean);



% now configure the synapse
x.set('*gmax',30);
x.reset;
x.Cell1.V = -61; % a little push to get it out of the synchronous state

x.t_end = 20e3;
V = x.integrate;
V(1:1e4,:) = [];

m(2) = xtools.V2metrics(V(:,2),'sampling_rate',10,'spike_threshold',-20);
m(1) = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);

COST = COST + B*xfit.binCost([500 2e3],m(1).burst_period) + B*xfit.binCost([500 2e3],m(2).burst_period);


% duty cycle
COST = COST + C*xfit.binCost([.1 .5],m(1).duty_cycle_mean) + C*xfit.binCost([.1 .5],m(2).duty_cycle_mean);



% phase

% find the burst starts and ends in the master cell 
try
	spike_times = xtools.findNSpikeTimes(V(:,1),xtools.findNSpikes(V(:,1)));
	spike_times2 = xtools.findNSpikeTimes(V(:,2),xtools.findNSpikes(V(:,2)));
	ibi = 3000;
	burst_ends = spike_times(find(diff(spike_times) > ibi));

	burst_starts = find(diff(spike_times) > ibi) + 1;
	burst_starts(burst_starts>length(spike_times)) = [];
	burst_starts = spike_times(burst_starts);

	if burst_ends(1) < burst_starts(1)
		burst_ends(1) = [];
	end

	n_ok_spikes = 0;
	n_wrong_spikes = 0;

	M = min([length(burst_starts) length(burst_ends)]);

	for i = 2:M-1
		n_wrong_spikes = n_wrong_spikes + sum(spike_times2 > burst_starts(i) &  spike_times2 < burst_ends(i));

		n_ok_spikes = n_ok_spikes + sum(spike_times2 > burst_ends(i-1) &  spike_times2 < burst_starts(i));

	end


	COST = COST + F*xfit.binCost([0, .1],n_wrong_spikes/(n_wrong_spikes + n_ok_spikes)) ...
   +F*xfit.binCost([.9, 1],n_ok_spikes/(n_wrong_spikes + n_ok_spikes));
catch
	COST = COST + 2*F;
end







if nargout > 1

	% measure behaviour with the synapse
	% and out of phase
	x.set('*gmax',30);
	x.reset;
	x.Cell2.V = -61;
	x.t_end = 20e3;
	[V, Ca] = x.integrate;
	Ca(1:1e4,:) = [];
	V(1:1e4,:) = [];
	metrics = xtools.V2metrics(V(:,1),'sampling_rate',10);
	metrics.Ca_mean = mean(Ca(:,1));
	metrics.Ca_std = std(Ca(:,1));
	metrics_new = xtools.V2metrics(V(:,2),'sampling_rate',10);
	metrics_new.Ca_mean = mean(Ca(:,2));
	metrics_new.Ca_std = std(Ca(:,2));
	% concat
	fn = fieldnames(metrics);
	for i = 1:length(fn)
		metrics.(fn{i}) = [metrics.(fn{i}) metrics_new.(fn{i})];
	end

	% measure behaviour with synapse
	% and starting in phase
	x.set('*gmax',30);
	x.reset;
	x.t_end = 20e3;
	[V, Ca] = x.integrate;
	Ca(1:1e4,:) = [];
	V(1:1e4,:) = [];
	metrics_new = xtools.V2metrics(V(:,1),'sampling_rate',10);
	metrics_new.Ca_mean = mean(Ca(:,1));
	metrics_new.Ca_std = std(Ca(:,1));

	for i = 1:length(fn)
		metrics.(fn{i}) = [metrics.(fn{i}) metrics_new.(fn{i})];
	end

	metrics_new = xtools.V2metrics(V(:,2),'sampling_rate',10);
	metrics_new.Ca_mean = mean(Ca(:,2));
	metrics_new.Ca_std = std(Ca(:,2));

	for i = 1:length(fn)
		metrics.(fn{i}) = [metrics.(fn{i}) metrics_new.(fn{i})];
	end




	% measure the calcium without the synapse
	x.set('*gmax',0);
	x.reset;
	x.t_end = 20e3;
	[V, Ca] = x.integrate;
	Ca(1:1e4,:) = [];
	V(1:1e4,:) = [];
	metrics_new = xtools.V2metrics(V(:,1),'sampling_rate',10);
	metrics_new.Ca_mean = mean(Ca(:,1));
	metrics_new.Ca_std = std(Ca(:,1));

	for i = 1:length(fn)
		metrics.(fn{i}) = [metrics.(fn{i}) metrics_new.(fn{i})];
	end

	metrics_new = xtools.V2metrics(V(:,2),'sampling_rate',10);
	metrics_new.Ca_mean = mean(Ca(:,2));
	metrics_new.Ca_std = std(Ca(:,2));

	for i = 1:length(fn)
		metrics.(fn{i}) = [metrics.(fn{i}) metrics_new.(fn{i})];
	end


end





