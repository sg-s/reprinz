function C = conditional_burster_cf(x,~,~)

x.reset;

C = 1e4;

x.t_end = 10e3;
x.dt = .1;
x.sim_dt = .1;

% turn synapse off
x.synapses.gbar = 0;

x.closed_loop = true;

x.integrate;
[V,Ca] = x.integrate;

% does the other neuron spike? 
m(2) = xtools.V2metrics(V(:,2),'sampling_rate',10,'spike_threshold',-20);
m(1) = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);

if m(2).firing_rate == 0
	return
else
	C = C - 1e3;
end

% firing rate should be higher than twice burst period
if m(2).firing_rate < 2/(m(1).burst_period*1e-3)
	return
else
	C = C - 1e3;
end

% should be firing regularly 
if m(2).isi_std/m(2).isi_mean > .01
	return
else
	C = C - 1e3;
end

% OK, it's OK w/o synapses
% now configure the synapse
C = 1e3;

x.synapses.gbar = 30;
x.reset;
x.integrate;
V = x.integrate;

m(2) = xtools.V2metrics(V(:,2),'sampling_rate',10,'spike_threshold',-20);
m(1) = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);

% follower cell should be bursting, not silent or spiking
if isnan(m(2).burst_period)
	return
elseif  m(2).isi_std/m(2).isi_mean < 1
	% probably spiking
	return
elseif m(2).firing_rate < m(1).firing_rate/2
	% not firing enough
	return
end

% find delays
Z = zscore(Ca(:,1:2));
d = finddelay(Z(:,1),Z(:,2));


% burst periods should be the same
% and the follower cell should burst regularly

C = 10*abs(m(2).burst_period/m(1).burst_period - 1) ...
    + m(2).ibi_std/m(2).ibi_mean ...
    + m(2).burst_period_std/m(2).burst_period ...
    + 20*abs(abs(d/(m(1).burst_period/2)) - 1);
