% this cost function can be used to find 
% neurons that tonically spike without any inhibitory
% input, and fire bursts between bursts of inhib.
% synaptic input
% 
function [C, Ca_burster, Ca_CB0, Ca_CB1] = conditional_burster_cf(x,~,~)

x.reset;

Ca_burster = NaN;
Ca_CB0 = NaN;
Ca_CB1 = NaN;

% don't change these!

x.dt = .1;
x.sim_dt = .1;

% turn synapse off
x.set('*gmax',0)


x.closed_loop = false;

x.t_end = 20e3;
V = x.integrate;
V(1:1e4,:) = [];

% does the other neuron spike? 
m(2) = xtools.V2metrics(V(:,2),'sampling_rate',10,'spike_threshold',-20);
m(1) = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);

if m(2).firing_rate == 0
	C = 1e4;
	return
end

C = 10*bin_cost([m(1).firing_rate/2, m(1).firing_rate*2],m(2).firing_rate);


if nargout == 0
	disp('Firing rate w/o synapses OK')
end

C = C+10*bin_cost([0,1e-2],m(2).isi_std/m(2).isi_mean);



% OK, it's OK w/o synapses
% now configure the synapse

x.synapses.gbar = 30;
x.reset;

x.t_end = 20e3;
V = x.integrate;
V(1:1e4,:) = [];

m(2) = xtools.V2metrics(V(:,2),'sampling_rate',10,'spike_threshold',-20);
m(1) = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);


C = C+10*bin_cost([m(1).firing_rate*.5, m(1).firing_rate*2],m(2).firing_rate);





% find the burst starts and ends in the master cell 
spike_times = xolotl.findNSpikeTimes(V(:,1),xolotl.findNSpikes(V(:,1)));
spike_times2 = xolotl.findNSpikeTimes(V(:,2),xolotl.findNSpikes(V(:,2)));
ibi = 3000;
burst_ends = spike_times(find(diff(spike_times) > ibi));

burst_starts = find(diff(spike_times) > ibi) + 1;
burst_starts(burst_starts>length(spike_times)) = [];
burst_starts = spike_times(burst_starts);

% time = (1:length(V))*1e-3*x.dt;
% plot(time,V(:,1)); hold on
% plot(time(burst_starts),V(burst_starts),'ro')
% plot(time(burst_ends),V(burst_ends),'ko')

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

C = C+10*bin_cost([0, .1],n_wrong_spikes/(n_wrong_spikes + n_ok_spikes)) ...
   +10*bin_cost([.9, 1],n_ok_spikes/(n_wrong_spikes + n_ok_spikes));

if nargout > 1

	% measure the calcium with the synapse
	x.synapses.gbar = 30;
	x.reset;
	x.t_end = 20e3;
	[~, Ca] = x.integrate;
	Ca(1:1e4,:) = [];
	Ca_burster = mean(Ca(:,1));
	Ca_CB1 = mean(Ca(:,2));


	% measure the calcium without the synapse
	x.synapses.gbar = 0;
	x.reset;
	[~, Ca] = x.integrate;
	Ca(1:1e4,:) = [];
	Ca_CB0 = mean(Ca(:,2));

end


end


function c = bin_cost(allowed_range,actual_value)


	w = (allowed_range(2) - allowed_range(1))/2;
	m = (allowed_range(2) + allowed_range(1))/2;

	if actual_value < allowed_range(1)
		d = m - actual_value;
		c = (1- (w/d));
	elseif actual_value > allowed_range(2)
		d = actual_value - m;
		c = (1- (w/d));
	else
		% no cost
		c = 0;
	end

end