% cost function that will be used to tune 
% a mulit-compartment neuron model


function C = realistic_cost_func(x,~,~)

;;;;;;;;     ;;;    ;;;;;;;;     ;;;    ;;     ;;  ;;;;;;  
;;     ;;   ;; ;;   ;;     ;;   ;; ;;   ;;;   ;;; ;;    ;; 
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;  ;;;; ;;;; ;;       
;;;;;;;;  ;;     ;; ;;;;;;;;  ;;     ;; ;; ;;; ;;  ;;;;;;  
;;        ;;;;;;;;; ;;   ;;   ;;;;;;;;; ;;     ;;       ;; 
;;        ;;     ;; ;;    ;;  ;;     ;; ;;     ;; ;;    ;; 
;;        ;;     ;; ;;     ;; ;;     ;; ;;     ;;  ;;;;;;  


% first set up the parmeters
g = x.get(x.find('Soma*gbar'));
comp_names = x.find('compartment');
for i = 1:length(comp_names)
	x.(comp_names{i}).set(x.(comp_names{i}).find('*gbar'),g);
end


% turn off NaV close to the soma
x.Soma.NaV.gbar = 0;
start_axon = floor(length(comp_names)/2);

for i = 1:start_axon
	x.(comp_names{i}).NaV.gbar = 0;
end


C = 1e9;
cost_vector = [1e4 1e3 3e2];

isi_range = [.02 .08];
ibi_range = [.5 1.5 ]; % seconds

ibi_cv_range = [0 .1]; % CV

CV_Ca_peak_period_range = [0 .1];
n_spikes_per_burst_range = [4 10];
max_min_V_within_burst = -40; 
cycle_period_range = [1.23 1.788]; % mean +/- 1 SD
burst_duration_range = [0.4490  0.7150]; % PD, sec
duty_cycle_range = [.3450 .4250];

% integrate, etc. 
[V,Ca] = x.integrate;
x.t_end = 15e3;
cutoff = floor(5e3/x.dt);
V_soma = V(cutoff:end,end);
V = V(cutoff:end,end-1);
T = length(V)*x.dt*1e-3; % seconds


;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;             ;;   
;;       ;;       ;;     ;; ;;       ;;           ;;;;   
;;       ;;       ;;     ;; ;;       ;;             ;;   
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;             ;;   
;;       ;;        ;;   ;;  ;;       ;;             ;;   
;;       ;;         ;; ;;   ;;       ;;             ;;   
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;     ;;;;;; 

% level 1 checks whether the cell is spiking or not
level_cost = cost_vector(1);

spike_times = psychopomp.findNSpikes(V,1e3);
n_spikes = sum(~isnan(spike_times));
if n_spikes > 0
	% penalize based on expected # of spikes
	min_expected_cycles = ((T)/cycle_period_range(2))-1;
	min_expected_spikes = min_expected_cycles*n_spikes_per_burst_range(1);
	if n_spikes < min_expected_spikes
		cost_vector(1) = level_cost*(n_spikes/min_expected_spikes);
	else
		% OK, zero cost
		cost_vector(1) = 0;
	end
else
	% give it a cost that grows with how hyperpolarized it is
	if max(V) > 0
		% always above zero
		cost_vector(1) = level_cost;
	else
		cost_vector(1) = (0 - min(V))*level_cost;
	end
	
end

spike_times = spike_times*x.dt*1e-3;

C = sum(cost_vector);

time = (1:length(V))*x.dt*1e-3;


if nargout == 0

	figure('outerposition',[300 300 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
	subplot(2,1,1); hold on
	plot(time,V,'k')
	ylabel('V_{terminal} (mV)')
	set(gca,'YLim',[-80 50])
	subplot(2,1,2); hold on
	plot(time,V_soma,'k')
	ylabel('V_{soma} (mV)')
	set(gca,'YLim',[-80 50])
	prettyFig('plw',1.5,'lw',1.5);
	

end

if cost_vector(1) > 0
	return
end

;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;           ;;;;;;;  
;;       ;;       ;;     ;; ;;       ;;          ;;     ;; 
;;       ;;       ;;     ;; ;;       ;;                 ;; 
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;           ;;;;;;;  
;;       ;;        ;;   ;;  ;;       ;;          ;;        
;;       ;;         ;; ;;   ;;       ;;          ;;        
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;    ;;;;;;;;; 

% OK, cell is spiking. now is it bursting, regularly? 

level_cost = cost_vector(2);

isis = diff(nonnans(spike_times));

% this way, every ISIs that is outside the range
% is considered a burst
ibis = isis(isis>isi_range(2));

% do we have the right # of IBIs?

n_ibi_cost =  .5*level_cost*bin_cost(sort(T./cycle_period_range),length(ibis));

if length(ibis) < 4
	cost_vector(2) = n_ibi_cost + .5*level_cost;
	C = sum(cost_vector);
	return
end


ibi_cv_cost = .5*level_cost*bin_cost(ibi_cv_range,cv(ibis));

cost_vector(2) = n_ibi_cost + ibi_cv_cost;
C = sum(cost_vector);

;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;           ;;;;;;;  
;;       ;;       ;;     ;; ;;       ;;          ;;     ;; 
;;       ;;       ;;     ;; ;;       ;;                 ;; 
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;           ;;;;;;;  
;;       ;;        ;;   ;;  ;;       ;;                 ;; 
;;       ;;         ;; ;;   ;;       ;;          ;;     ;; 
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;     ;;;;;;;  

level_cost = cost_vector(3);

% now look at duty cycle, etc. 

% find all burst starts and stops
burst_starts = spike_times(find(isis > ibi_range(1))+1);

if ~isempty(burst_starts)

	burst_ends = NaN*burst_starts;
	spike_times = nonnans(spike_times);
	% find the last spike after each burst start
	for i = 2:length(burst_starts)-1

		temp = max(spike_times(spike_times > burst_starts(i) & (spike_times < burst_starts(i+1))));
		if ~isempty(temp)
			burst_ends(i) = temp;
		end

	end

	burst_ends = burst_ends(2:end-1);
	burst_starts = burst_starts(2:end-1);


	% measure # of spikes/burst
	nspb = 0*burst_starts;
	for i = 1:length(burst_starts)
		nspb(i) = sum(spike_times > burst_starts(i) & spike_times < burst_ends(i));
	end
	if cv(nspb) < .1
		nspb = mean(nspb);

		% deliberately over weight this cost
		nspb_cost = (level_cost)*bin_cost(n_spikes_per_burst_range,nspb);
	else
		% too variable, give it a high cost
		nspb_cost = level_cost;
	end


	% measure burst frequencies of neuron
	all_periods = mean(diff(burst_starts));
	all_durations = nanmean(burst_ends - burst_starts);

	period_cost = (level_cost/4)*bin_cost(cycle_period_range,all_periods);
	duration_cost = (level_cost/4)*bin_cost(burst_duration_range,nanmean(all_durations));
	all_duty_cycle = nanmean(all_durations)/all_periods;
	dc_cost = (level_cost/4)*bin_cost(duty_cycle_range,all_duty_cycle);


	% now also check for the real ISIs
	true_isis = isis(isis<ibi_range(1));
	isi_cost = 0*true_isis;
	for i = 1:length(true_isis)
		isi_cost(i) = bin_cost(isi_range,true_isis(i));
	end
	isi_cost = sum(isi_cost)*(level_cost/4);


	if ~nargout
		disp(['Period cost = ' oval(period_cost)])
		disp(['Duration cost = ' oval(duration_cost)])
		disp(['Duty cycle cost = ' oval(dc_cost)])
		disp(['#spikes/burst  cost = ' oval(nspb_cost)])
		disp(['ISI  cost = ' oval(isi_cost)])
	end

	cost_vector(3) =  period_cost + duration_cost + dc_cost + nspb_cost + isi_cost;


end

	C = sum(cost_vector(:));

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

end