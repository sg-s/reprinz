% this cost function can be used to find 
% neurons that tonically spike without any inhibitory
% input, and fire bursts between bursts of inhib.
% synaptic input
% 
function [C, Ca_with, Ca_without] = hco_cost_func(x,~,~)

x.reset;

% mirror gbars in Cell2
x.set('Cell2*gbar',x.get('Cell1*gbar'));

Ca_burster = NaN;
Ca_CB0 = NaN;
Ca_CB1 = NaN;

% don't change these!

x.dt = .1;
x.sim_dt = .1;

% turn synapse off. now, we want both neurons to be either silent
% or tonically spiking
x.set('*gmax',0)

x.closed_loop = false;

x.t_end = 20e3;
V = x.integrate;
V(1:1e4,:) = [];

m(2) = xtools.V2metrics(V(:,2),'sampling_rate',10,'spike_threshold',-20);
m(1) = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);



if m(2).firing_rate == 0 & m(1).firing_rate == 0
	% both silent, nice
elseif m(1).isi_std/m(1).isi_mean > .01
	% bursting 
	C = 1e5;
	return
else
	%disp('Spiking?')
	% spiking is OK
end


% OK, it's OK w/o synapses
% now configure the synapse

x.set('*gmax',30);
x.reset;
x.Cell1.V = 0;

x.t_end = 20e3;
V = x.integrate;
V(1:1e4,:) = [];

m(2) = xtools.V2metrics(V(:,2),'sampling_rate',10,'spike_threshold',-20);
m(1) = xtools.V2metrics(V(:,1),'sampling_rate',10,'spike_threshold',-20);

% make sure both cells are spiking 
if m(1).firing_rate == 0 | m(2).firing_rate == 0
	C = 1e4;
	%disp('Not spiking enough, aborting...')
	return
end



% make sure both cells have similar firing rates
df = abs(m(1).firing_rate  - m(2).firing_rate)/(m(2).firing_rate+m(1).firing_rate);
df = df/2;

C = bin_cost([0,1e-1],df);
if C > 0
	C = C*1e3;
	%disp('Firing rates of two neurons mismatched...')
	return
end


% make sure that the cells are bursting
% force regular bursting
C = 0;
if m(1).n_spikes_per_burst_std/m(1).n_spikes_per_burst_mean > 1e-2 || m(2).n_spikes_per_burst_std/m(2).n_spikes_per_burst_mean > 1e-2
	C = C + 100;
	%disp('Irregular burster...')
end

if isnan(m(1).burst_period)
	C = C + 100;
end

if isnan(m(2).burst_period)
	C = C + 100;
end

% make sure they're not 1 spike bursters
C = C + 100*bin_cost([4 30],m(1).n_spikes_per_burst_mean);
C = C  + 100*bin_cost([4 30],m(2).n_spikes_per_burst_mean);
if m(1).n_spikes_per_burst_mean < 2 | m(2).n_spikes_per_burst_mean < 2
	C = C*10; % weight this a little extra
	%disp('Single-spike bursters...')
end

% check that the burst periods are approx eq.
df = abs(m(1).burst_period  - m(2).burst_period)/(m(2).burst_period+m(1).burst_period);
df = df/2;

C = C+20*bin_cost([0,1e-1],df);

if C > 0
	%disp('bursting not OK')
	return
end

% find the burst starts and ends in the master cell 
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

C = C+10*bin_cost([0, .1],n_wrong_spikes/(n_wrong_spikes + n_ok_spikes)) ...
   +10*bin_cost([.9, 1],n_ok_spikes/(n_wrong_spikes + n_ok_spikes));

if nargout > 1

	% measure the calcium with the synapse
	x.set('*gmax',30);
	x.reset;
	x.Cell2.V = 0;
	x.t_end = 20e3;
	[~, Ca] = x.integrate;
	Ca(1:1e4,:) = [];
	Ca_with = [mean(Ca(:,1)); mean(Ca(:,2))];

	% measure the calcium without the synapse
	x.set('*gmax',0);
	x.reset;
	[~, Ca] = x.integrate;
	Ca(1:1e4,:) = [];
	Ca_without = [mean(Ca(:,1)); mean(Ca(:,2))];

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