% this cost function is meant to find networks
% that are within the experimental range 
% that is observed for STG networks
% in a number of parameters as defined in Prinz et al 2004

function [C, cost_vector] = STG_cost_function(x)

% first, make the default cost_vector
cost_vector = zeros(4,3);
cost_vector(1,:) = 1e5;
cost_vector(2,:) = 1e4;
cost_vector(3,:) = 1e3;
cost_vector(4,:) = 1e2;
C = sum(cost_vector(:));

n_spikes_per_burst_range = [3 30; 3 30; 3 30];

% first, we simulate the network 
[V,Ca] = x.integrate;

% throw away 10 seconds of transient 
cutoff = floor(10e3/x.dt);
V = V(cutoff:end,:);
Ca = Ca(cutoff:end,:);

time = (1:length(V))*x.dt;

if nargout == 0
	figure('outerposition',[300 300 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
	for i = 3:-1:1
		ax(i) = subplot(3,1,i); hold on
		plot(time,V(:,i),'k')
	end
end



cycle_period_range = [1.23 1.788]; % mean +/- 1 SD
PD_burst_duration_range = [0.4490  0.7150]; % sec
LP_burst_duration_range = [0.2860 .5120]; % sec
PY_burst_duration_range = [0.38 .68]; % sec

% define gaps A_end -> B_start
Gap_PD_LP_range = [.112 .33];
Gap_LP_PY_range = [-.121 -.001];

% delays are X_start -> B_start
Delay_PD_LP_range = [.6340 .9720]; % seconds
Delay_PD_PY_range = [.9250 1.3570]; % seconds


% duty cycles
PD_duty_cycle_range = [.3450 .4250];
LP_duty_cycle_range = [.2050 .3230];
PY_duty_cycle_range = [.2940 .4020];


;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;             ;;   
;;       ;;       ;;     ;; ;;       ;;           ;;;;   
;;       ;;       ;;     ;; ;;       ;;             ;;   
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;             ;;   
;;       ;;        ;;   ;;  ;;       ;;             ;;   
;;       ;;         ;; ;;   ;;       ;;             ;;   
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;     ;;;;;; 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% check that every cell is not silent
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for i = 3:-1:1
	spike_times(:,i) = psychopomp.findNSpikes(V(:,i),1e3);
	if sum(~isnan(spike_times(:,i))) > 10
		cost_vector(1,i) = 0;
	end
end

C = sum(cost_vector(:));
if C >= 3e5;
	return
end

;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;           ;;;;;;;  
;;       ;;       ;;     ;; ;;       ;;          ;;     ;; 
;;       ;;       ;;     ;; ;;       ;;                 ;; 
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;           ;;;;;;;  
;;       ;;        ;;   ;;  ;;       ;;          ;;        
;;       ;;         ;; ;;   ;;       ;;          ;;        
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;    ;;;;;;;;; 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% cell is bursting, with an OK # of spikes/burst
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

spike_times = spike_times*x.dt*1e-3;


% now we check to see if the minimum calcium is less 
% that 10% of the peak 
not_bursting = false;
for i = 1:3
	if min(Ca(:,i)) < .1*max(Ca(:,i))
		cost_vector(2,i) = cost_vector(2,i) - .5e4;
	else
		not_bursting = true;
	end
end

% early exit
if not_bursting
	C = sum(cost_vector(:));
	return
end


% now look at the # of spike/burst in each cell
% for that we first have to get the spiketimes/burst
Ca_peak_times = NaN(100,3);
Ca_min_times = NaN(100,3);
Ca_peak_values = NaN(100,3);

for i = 1:3
	[ons,offs] = computeOnsOffs(Ca(:,i) > mean(Ca(:,i))/2);
	% find peaks b/w ons and offs
	for j = 1:length(ons)
		[Ca_peak_values(j,i),idx] = max(Ca(ons(j):offs(j),i));
		Ca_peak_times(j,i) = ons(j) + idx - 1;
	end
	% find mins b/w offs and next ons 
	for j = 1:length(ons)-1
		[~,idx] = max(-Ca(offs(j):ons(j+1),i));
		Ca_min_times(j,i) = offs(j) + idx - 1;
	end
end


% convert these things into seconds
Ca_peak_times = Ca_peak_times*x.dt*1e-3;
Ca_min_times = Ca_min_times*x.dt*1e-3;

for i = 1:3
	n_bursts = sum(~isnan(Ca_min_times(:,i)));
	frac_cost = cost_vector(2,i)/(n_bursts-1);
	for j = 1:n_bursts-1
		spikes_in_this_burst = spike_times(Ca_min_times(j,i) <= spike_times(:,i) & Ca_min_times(j+1,i) >= spike_times(:,i),i);
		if (length(spikes_in_this_burst) > n_spikes_per_burst_range(i,1)) & (length(spikes_in_this_burst) <n_spikes_per_burst_range(i,2))
			% in range, relax cost
			cost_vector(2,i) = cost_vector(2,i) - frac_cost;
		else
			% disp(['not in range: ' mat2str(length(spikes_in_this_burst))])
		end
	end
end


C = sum(cost_vector(:));


if C >= 3e4
	return
end

;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;           ;;;;;;;  
;;       ;;       ;;     ;; ;;       ;;          ;;     ;; 
;;       ;;       ;;     ;; ;;       ;;                 ;; 
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;           ;;;;;;;  
;;       ;;        ;;   ;;  ;;       ;;                 ;; 
;;       ;;         ;; ;;   ;;       ;;          ;;     ;; 
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;     ;;;;;;;  

level_cost = 1e3;

% check for depolarization block
for i = 1:3
	n_bursts = sum(~isnan(Ca_min_times(:,i)));
	frac_cost = level_cost/((n_bursts-1));
	for j = 1:n_bursts-1
		spikes_in_this_burst = spike_times(Ca_min_times(j,i) <= spike_times(:,i) & Ca_min_times(j+1,i) >= spike_times(:,i),i);
		

		isis_in_burst = diff(spikes_in_this_burst);

		if max(isis_in_burst)/mean(isis_in_burst) < 3
			cost_vector(3,i) = cost_vector(3,i) - frac_cost;
		end

	end
end

C = sum(cost_vector(:));

if C >= 3e3
	return
end

;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;          ;;        
;;       ;;       ;;     ;; ;;       ;;          ;;    ;;  
;;       ;;       ;;     ;; ;;       ;;          ;;    ;;  
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;          ;;    ;;  
;;       ;;        ;;   ;;  ;;       ;;          ;;;;;;;;; 
;;       ;;         ;; ;;   ;;       ;;                ;;  
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;          ;;  

% here we check burst frequencies and duty cycles and make sure that they are all within the experimental range

level_cost = 1e2; % we're measuring ten things, so total level cost is 1e2

% measure burst frequencies of all neurons
for i = 3:-1:1
	all_freq(i) = mean(diff(nonnans(Ca_min_times(:,i))));
end

% penalize them based on how far they are from the 
% experimental range 
w = (cycle_period_range(2) - cycle_period_range(1))/2;
m = (cycle_period_range(2) + cycle_period_range(1))/2;
for i = 1:3
	if all_freq(i) < cycle_period_range(1)
		d = m - all_freq(i);
		cost_vector(4,i) = (1- (w/d))*level_cost;
	elseif all_freq(i) > cycle_period_range(2)
		d = all_freq(i) - m;
		cost_vector(4,i) = (1- (w/d))*level_cost;
	else
		% no cost
		cost_vector(4,i) = 0;
	end
end

% add a penalty if the frequencies don't match up
if max(all_freq)/min(all_freq) > 1.01
	cost_vector(4,:) = cost_vector(4,:) + level_cost;
end


C = sum(cost_vector(:));
if C >= 3e2
	return
end


return


% first find burst metrics for all three neurons
for i = 3:-1:1
	metrics(i) = psychopomp.findBurstMetrics(V(:,i),Ca(:,i),Inf, Inf,0);

	metrics(i).spike_times = metrics(i).spike_times*x.dt;
	metrics(i).Ca_peaks = metrics(i).Ca_peaks*x.dt;
	metrics(i).Ca_mins = metrics(i).Ca_mins*x.dt;
	metrics(i).last_spike_loc = metrics(i).last_spike_loc*x.dt;
	metrics(i).first_spike_loc = metrics(i).first_spike_loc*x.dt;


	% early exit -- check if everything is OK
	if metrics(i).burst_metrics(10) ~= 0
		C = fail_cost;
		return
	end

	if isempty(metrics(i).Ca_mins) || isempty(metrics(i).first_spike_loc) || isempty(metrics(i).last_spike_loc) 
		C = fail_cost;
		return
	end

end




% now measure the neuron-specific things

try

	for i = 1:3
		% measure the cost of the cycle period for this neuron
		C = C +  basin_cost(cycle_period_range, metrics(i).burst_metrics(1)*x.dt*1e-3);

		% measure the duration of this neuron
		this_duration = [];
		for j = (length(metrics(i).Ca_mins)-1):-1:1
			if isnan(metrics(i).Ca_mins(j))
				continue
			end
			if isnan(metrics(i).Ca_mins(j+1))
				continue
			end
			these_spikes = metrics(i).spike_times(metrics(i).spike_times > metrics(i).Ca_mins(j) & metrics(i).spike_times < metrics(i).Ca_mins(j+1));
			this_duration(j) = these_spikes(end) - these_spikes(1);
		end
		this_duration = 1e-3*mean(nonzeros(this_duration));


		if i == 1
			% PD neuron

			% measure PD neuron duration (1st spike -- last spike)
			C =  C + basin_cost(PD_burst_duration_range,this_duration);
		end

		if i == 2
			% LP neuron

			% measure LP neuron duration (1st spike -- last spike)
			C =  C + basin_cost(LP_burst_duration_range,this_duration);
		end

		if i == 3
			% PY neuron

			% measure PY neuron duration (1st spike -- last spike)
			C =  C + basin_cost(PY_burst_duration_range,this_duration);
		end

	end  % end neuron-specific loop

	PD_ends = metrics(1).last_spike_loc;
	rm_this = PD_ends == 0;
	PD_ends = metrics(1).Ca_peaks(:) + PD_ends(:);
	PD_ends(rm_this) = [];

	PD_starts = metrics(1).first_spike_loc;
	rm_this = PD_starts == 0;
	PD_starts = metrics(1).Ca_peaks(:) + PD_starts(:);
	PD_starts(rm_this) = [];


	LP_starts = metrics(2).first_spike_loc;
	rm_this = LP_starts == 0;
	LP_starts = metrics(2).Ca_peaks(:) + LP_starts(:);
	LP_starts(rm_this) = [];

	LP_ends = metrics(2).last_spike_loc;
	rm_this = LP_ends == 0;
	LP_ends = metrics(2).Ca_peaks(:) + LP_ends(:);
	LP_ends(rm_this) = [];


	PY_starts = metrics(3).first_spike_loc;
	rm_this = PY_starts == 0;
	PY_starts = metrics(3).Ca_peaks(:) + PY_starts(:);
	PY_starts(rm_this) = [];

	PY_ends = metrics(3).last_spike_loc;
	rm_this = PY_ends == 0;
	PY_ends = metrics(3).Ca_peaks(:) + PY_ends(:);
	PY_ends(rm_this) = [];

catch
	C = fail_cost;
	return
end


if nargout == 0
	% mark the burst ends and starts
	plot(ax(1),[PD_starts PD_starts],[-80 60],'g','LineWidth',3)
	plot(ax(1),[PD_ends PD_ends],[-80 60],'r','LineWidth',3)

	plot(ax(2),[LP_starts LP_starts],[-80 60],'g','LineWidth',3)
	plot(ax(2),[LP_ends LP_ends],[-80 60],'r','LineWidth',3)

	plot(ax(3),[PY_starts PY_starts],[-80 60],'g','LineWidth',3)
	plot(ax(3),[PY_ends PY_ends],[-80 60],'r','LineWidth',3)
end

try

% measure gaps from PD to LP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_gaps = NaN*PD_ends;

for i = 2:length(all_gaps)-1
	[~,idx] = min(abs(PD_ends(i) - LP_starts));
	all_gaps(i) = (LP_starts(idx) - PD_ends(i))*1e-3;
end

C = C + basin_cost(Gap_PD_LP_range, nanmean(all_gaps));

% measure gaps from LP to PY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_gaps = NaN*LP_ends;

for i = 2:length(all_gaps)-1
	[~,idx] = min(abs(LP_ends(i) - PY_starts));
	all_gaps(i) = (PY_starts(idx) - LP_ends(i))*1e-3;
end

C = C + basin_cost(Gap_LP_PY_range, nanmean(all_gaps));

;;;;;;;;  ;;;;;;;; ;;          ;;;    ;;    ;;  ;;;;;;  
;;     ;; ;;       ;;         ;; ;;    ;;  ;;  ;;    ;; 
;;     ;; ;;       ;;        ;;   ;;    ;;;;   ;;       
;;     ;; ;;;;;;   ;;       ;;     ;;    ;;     ;;;;;;  
;;     ;; ;;       ;;       ;;;;;;;;;    ;;          ;; 
;;     ;; ;;       ;;       ;;     ;;    ;;    ;;    ;; 
;;;;;;;;  ;;;;;;;; ;;;;;;;; ;;     ;;    ;;     ;;;;;;  

% measure delay from PD -> LP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_delays = NaN*PD_starts;

for i = 2:length(all_delays)-1
	[~,idx] = min(abs(PD_starts(i) - LP_starts));
	all_delays(i) = (LP_starts(idx) - PD_starts(i))*1e-3;
end

C = C + basin_cost(Delay_PD_LP_range, nanmean(all_delays));

% measure delay from PD -> PY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_delays = NaN*PD_starts;

for i = 2:length(all_delays)-1
	[~,idx] = min(abs(PD_starts(i) - PY_starts));
	all_delays(i) = (PY_starts(idx) - PD_starts(i))*1e-3;
end

C = C + basin_cost(Delay_PD_PY_range, nanmean(all_delays));


;;;;;;;;  ;;     ;; ;;;;;;;; ;;    ;; 
;;     ;; ;;     ;;    ;;     ;;  ;;  
;;     ;; ;;     ;;    ;;      ;;;;   
;;     ;; ;;     ;;    ;;       ;;    
;;     ;; ;;     ;;    ;;       ;;    
;;     ;; ;;     ;;    ;;       ;;    
;;;;;;;;   ;;;;;;;     ;;       ;;    

 ;;;;;;  ;;    ;;  ;;;;;;  ;;       ;;;;;;;; 
;;    ;;  ;;  ;;  ;;    ;; ;;       ;;       
;;         ;;;;   ;;       ;;       ;;       
;;          ;;    ;;       ;;       ;;;;;;   
;;          ;;    ;;       ;;       ;;       
;;    ;;    ;;    ;;    ;; ;;       ;;       
 ;;;;;;     ;;     ;;;;;;  ;;;;;;;; ;;;;;;;; 

C = C + basin_cost(PD_duty_cycle_range,metrics(1).burst_metrics(9));
C = C + basin_cost(LP_duty_cycle_range,metrics(2).burst_metrics(9));
C = C + basin_cost(PY_duty_cycle_range,metrics(3).burst_metrics(9));


% penalize depolarization block in LP

for i = 1:length(LP_starts)
	these_spikes = metrics(2).spike_times(metrics(2).spike_times >= LP_starts(i) & (metrics(2).spike_times <= LP_ends(i)));

	allowed_isi_range = [0 3*(these_spikes(end) - these_spikes(1))/length(these_spikes)];

	C = C + basin_cost(allowed_isi_range,max(diff(these_spikes)));

end

% In PD
for i = 1:length(PD_starts)
	these_spikes = metrics(2).spike_times(metrics(2).spike_times >= PD_starts(i) & (metrics(2).spike_times <= PD_ends(i)));

	allowed_isi_range = [0 3*(these_spikes(end) - these_spikes(1))/length(these_spikes)];

	C = C + basin_cost(allowed_isi_range,max(diff(these_spikes)));

end




catch

	C = fail_cost;
	return
end



function c = basin_cost(allowed_range, actual_value)

	c = 0;
	if actual_value < allowed_range(1)
		c = (allowed_range(2) - actual_value)^2/(allowed_range(2) - allowed_range(1));
		return
	elseif  actual_value > allowed_range(2)
		c = (allowed_range(2) - actual_value)^2/(allowed_range(2) - allowed_range(1));
		return
	else
		return 
	end


end % end basin_cost

end % end function 