% this cost function is meant to find networks
% that are within the experimental range 
% that is observed for STG networks
% in a number of parameters as defined in Prinz et al 2004

function [C, cost_vector] = STG_cost_function(x)

% first, make the default cost_vector
cost_vector = zeros(6,3);
cost_vector(1,:) = 1e5; % is it bursting?  
cost_vector(2,:) = 1e4;
cost_vector(3,:) = 1e3;
cost_vector(4,:) = 1e2; % burst frequency
cost_vector(5,:) = 1e2; % duty cycle
cost_vector(6,:) = 1e2; % duration 


C = sum(cost_vector(:));


;;;;;;;;     ;;;    ;;;;;;;;     ;;;    ;;     ;;  ;;;;;;  
;;     ;;   ;; ;;   ;;     ;;   ;; ;;   ;;;   ;;; ;;    ;; 
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;  ;;;; ;;;; ;;       
;;;;;;;;  ;;     ;; ;;;;;;;;  ;;     ;; ;; ;;; ;;  ;;;;;;  
;;        ;;;;;;;;; ;;   ;;   ;;;;;;;;; ;;     ;;       ;; 
;;        ;;     ;; ;;    ;;  ;;     ;; ;;     ;; ;;    ;; 
;;        ;;     ;; ;;     ;; ;;     ;; ;;     ;;  ;;;;;;  


CV_Ca_peak_period_range = [0 .1];

n_spikes_per_burst_range = NaN(2,3);
n_spikes_per_burst_range(:,1) = [4 30];
n_spikes_per_burst_range(:,2) = [4 30];
n_spikes_per_burst_range(:,3) = [4 30];

max_min_V_within_burst = -40; 

cycle_period_range = [1.23 1.788]; % mean +/- 1 SD

burst_duration_range = NaN(2,3);
burst_duration_range(:,1) = [0.4490  0.7150]; % PD, sec
burst_duration_range(:,2) = [0.2860 .5120]; % LP, sec
burst_duration_range(:,3) = [0.38 .68]; % PY, sec


% duty cycles
duty_cycle_range = NaN(2,3);
duty_cycle_range(:,1) = [.3450 .4250];
duty_cycle_range(:,2) = [.2050 .3230];
duty_cycle_range(:,3) = [.2940 .4020];



% define gaps A_end -> B_start
Gap_PD_LP_range = [.112 .33];
Gap_LP_PY_range = [-.121 -.001];

% delays are X_start -> B_start
Delay_PD_LP_range = [.6340 .9720]; % seconds
Delay_PD_PY_range = [.9250 1.3570]; % seconds

;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;             ;;   
;;       ;;       ;;     ;; ;;       ;;           ;;;;   
;;       ;;       ;;     ;; ;;       ;;             ;;   
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;             ;;   
;;       ;;        ;;   ;;  ;;       ;;             ;;   
;;       ;;         ;; ;;   ;;       ;;             ;;   
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;     ;;;;;; 

level_cost = 1e5;


[V,Ca] = x.integrate;
cutoff = floor(5e3/x.dt);
V = V(cutoff:end,:);
Ca = Ca(cutoff:end,:);


for i = 3:-1:1
	spike_times(:,i) = psychopomp.findNSpikes(V(:,i),1e3);
	n_spikes = sum(~isnan(spike_times(:,i)));
	if n_spikes > 0
		% penalize based on expected # of spikes
		min_expected_cycles = (x.t_end*1e-3)/cycle_period_range(1);
		min_expected_spikes = min_expected_cycles*n_spikes_per_burst_range(1,i);
		if n_spikes < min_expected_spikes
			cost_vector(1,i) = level_cost*(n_spikes/min_expected_spikes);
		else
			% OK, zero cost
			cost_vector(1,i) = 0;
		end
	else
		% give it a cost that grows with how hyperpolarized it is
		if max(V(:,i)) > 0
			% always above zero
			cost_vector(1,i) = level_cost;
		else
			cost_vector(1,i) = (0 - min(V(:,i)))*level_cost;
		end
		
	end
end
spike_times = spike_times*x.dt*1e-3;


time = (1:length(V))*x.dt*1e-3;

if nargout == 0
	figure('outerposition',[300 300 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
	for i = 3:-1:1
		ax(i) = subplot(3,1,i); hold on
		plot(time,V(:,i),'k')
	end
end

C = sum(cost_vector(:));
if any(cost_vector(1,:) > 1)
	return
end


;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;           ;;;;;;;  
;;       ;;       ;;     ;; ;;       ;;          ;;     ;; 
;;       ;;       ;;     ;; ;;       ;;                 ;; 
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;           ;;;;;;;  
;;       ;;        ;;   ;;  ;;       ;;          ;;        
;;       ;;         ;; ;;   ;;       ;;          ;;        
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;    ;;;;;;;;; 

level_cost = 1e4;
cost_vector(2,i) = 0;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% cell is bursting, with an OK # of spikes/burst
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% step 1 of 3. now we check to see if the minimum calcium is less 
% that 10% of the peak 

for i = 1:3
	cost_vector(2,i) = bin_cost([0,.1],min(Ca(:,i))/max(Ca(:,i)))*level_cost/3;
end

% early exit
if any(cost_vector(2,:) > 1)
	C = sum(cost_vector(:));
	return
end


% figure out burst starts and stops based on calcium
% signals 
Ca_peak_times = NaN(100,3);
Ca_min_times = NaN(100,3);
Ca_peak_values = NaN(100,3);


for i = 1:3
	[ons,offs] = computeOnsOffs(Ca(:,i) > mean(Ca(:,i))/2);

	if length(ons) ~= length(offs)
		return
	end

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


if nargout == 0

	figure('outerposition',[300 300 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
	for i = 3:-1:1
		ax(i) = subplot(3,1,i); hold on
		plot(time,Ca(:,i),'k')
		temp = nonnans(Ca_peak_times(:,i));
		plot(time(temp),Ca(temp,i),'ro')

		temp = nonnans(Ca_min_times(:,i));
		plot(time(temp),Ca(temp,i),'bo')
	end

end

% convert these things into seconds
Ca_peak_times = Ca_peak_times*x.dt*1e-3;
Ca_min_times = Ca_min_times*x.dt*1e-3;

% step 2 of 3 -- check that Ca peaks occur regularly

for i = 1:3
	this_cv = cv(diff(nonnans(Ca_peak_times(:,i))));
	cost_vector(2,i) = cost_vector(2,i)+  level_cost*bin_cost(CV_Ca_peak_period_range,this_cv)/3;
end


% let's also store all the durations 
% and the # of spikes per burst
all_durations = NaN(100,3);
n_spikes_in_burst = NaN(100,3);
all_burst_starts = NaN(100,3);
all_burst_ends = NaN(100,3);

% check for depolarization block
for i = 1:3
	n_bursts = sum(~isnan(Ca_min_times(:,i)));
	for j = 1:n_bursts-1

		a = Ca_min_times(j,i);
		z = (Ca_min_times(j+1,i) - a)*.9 + a;

		spikes_in_this_burst = spike_times(a <= spike_times(:,i) & z >= spike_times(:,i),i);

		n_spikes_in_burst(j,i) = length(spikes_in_this_burst);

		% measure the duration of this burst
		if length(spikes_in_this_burst) > 1
			all_durations(j,i) = spikes_in_this_burst(end) - spikes_in_this_burst(1);
			all_burst_starts(j,i) = spikes_in_this_burst(1);
			all_burst_ends(j,i) = spikes_in_this_burst(end);
		else
			all_durations(j,i) = 0;
			continue
		end
	end
end

% penalize things that are outside the allowed
% # of spikes/burst


for i = 1:3
	this_n_spikes = nonnans(n_spikes_in_burst(:,i));
	frac_cost = (level_cost/3)/length(this_n_spikes);
	for j = 1:length(this_n_spikes)
		if (this_n_spikes(j) <= n_spikes_per_burst_range(2,i)) && (this_n_spikes(j) >= n_spikes_per_burst_range(1,i))
		else
			cost_vector(2,i) = cost_vector(2,i) + frac_cost;
		end
	end
end


C = sum(cost_vector(:));


if any(cost_vector(2,:) > 1)
	return
end


;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;           ;;;;;;;  
;;       ;;       ;;     ;; ;;       ;;          ;;     ;; 
;;       ;;       ;;     ;; ;;       ;;                 ;; 
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;           ;;;;;;;  
;;       ;;        ;;   ;;  ;;       ;;                 ;; 
;;       ;;         ;; ;;   ;;       ;;          ;;     ;; 
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;     ;;;;;;;  

% level 3 is all about the depolarization block. we don't
% want that. we hate depol. blocks. bad block. 

level_cost = 1e3;
cost_vector(3,:) = 0;

% check 1/2: measure V inbetween spikes

for i = 1:3
	n_bursts = sum(~isnan(n_spikes_in_burst(:,i)));
	frac_cost = (level_cost/2)/n_bursts;
	for j = 1:n_bursts
		spikes_in_this_burst = spike_times(all_burst_starts(j,i) <= spike_times(:,i) & all_burst_ends(j,i) >= spike_times(:,i),i);
		between_spikes = spikes_in_this_burst(1:end-1) + diff(spikes_in_this_burst)/2;

		max_V = max(V(round(between_spikes*1e3/x.dt),i));
		cost_vector(3,i) = cost_vector(3,i) + bin_cost([-60, -35],max_V)*frac_cost;

	end
end


% check 2/2 ISIs during bursts
for i = 1:3
	n_bursts = sum(~isnan(n_spikes_in_burst(:,i)));
	frac_cost = (level_cost/2)/n_bursts;
	for j = 1:n_bursts
		%expected_mean_isi = all_durations(j,i)/(n_spikes_in_burst(j,i)-1);
		spikes_in_this_burst = spike_times(all_burst_starts(j,i) <= spike_times(:,i) & all_burst_ends(j,i) >= spike_times(:,i),i);
		isis_in_this_burst = diff(spikes_in_this_burst);



		isi_ratio = max(isis_in_this_burst)/mean(setdiff(isis_in_this_burst,max(isis_in_this_burst)));
		cost_vector(3,i) = cost_vector(3,i) + bin_cost([0 5],isi_ratio)*frac_cost;

	end
end


C = sum(cost_vector(:));



if any(cost_vector(3,:) > 1)
	return
end

;;       ;;;;;;;; ;;     ;; ;;;;;;;; ;;          ;;        
;;       ;;       ;;     ;; ;;       ;;          ;;    ;;  
;;       ;;       ;;     ;; ;;       ;;          ;;    ;;  
;;       ;;;;;;   ;;     ;; ;;;;;;   ;;          ;;    ;;  
;;       ;;        ;;   ;;  ;;       ;;          ;;;;;;;;; 
;;       ;;         ;; ;;   ;;       ;;                ;;  
;;;;;;;; ;;;;;;;;    ;;;    ;;;;;;;; ;;;;;;;;          ;;  

% here we check that neuron-specific parameters like the period, duration and duty cycle are within the experimental range. 

level_cost = 1e2; 
cost_vector(4,:) = 0;

% measure burst frequencies of all neurons
for i = 3:-1:1
	all_periods(i) = mean(diff(nonnans(Ca_min_times(:,i))));
end

% penalize them based on how far they are from the 
% experimental range 

for i = 1:3
	cost_vector(4,i) = cost_vector(4,i) + bin_cost(cycle_period_range,all_periods(i));
end

% add a penalty if the frequencies don't match up
sync_cost = bin_cost([.95 1.05],max(all_periods)/min(all_periods));
cost_vector(4,:) = cost_vector(4,:) + sync_cost;

% durations 
for i = 1:3
	cost_vector(4,i) = cost_vector(4,i) + level_cost*bin_cost(burst_duration_range(:,i),nanmean(all_durations(:,i)));
end


% duty cycle

all_duty_cycle = nanmean(all_durations)./all_periods;
for i = 1:3
	cost_vector(6,i) = level_cost*bin_cost(duty_cycle_range(:,i),all_duty_cycle(i));
end


C = sum(cost_vector(:));

return


% % first find burst metrics for all three neurons
% for i = 3:-1:1
% 	metrics(i) = psychopomp.findBurstMetrics(V(:,i),Ca(:,i),Inf, Inf,0);

% 	metrics(i).spike_times = metrics(i).spike_times*x.dt;
% 	metrics(i).Ca_peaks = metrics(i).Ca_peaks*x.dt;
% 	metrics(i).Ca_mins = metrics(i).Ca_mins*x.dt;
% 	metrics(i).last_spike_loc = metrics(i).last_spike_loc*x.dt;
% 	metrics(i).first_spike_loc = metrics(i).first_spike_loc*x.dt;


% 	% early exit -- check if everything is OK
% 	if metrics(i).burst_metrics(10) ~= 0
% 		C = fail_cost;
% 		return
% 	end

% 	if isempty(metrics(i).Ca_mins) || isempty(metrics(i).first_spike_loc) || isempty(metrics(i).last_spike_loc) 
% 		C = fail_cost;
% 		return
% 	end

% end




% % now measure the neuron-specific things

% try

% 	for i = 1:3
% 		% measure the cost of the cycle period for this neuron
% 		C = C +  basin_cost(cycle_period_range, metrics(i).burst_metrics(1)*x.dt*1e-3);

% 		% measure the duration of this neuron
% 		this_duration = [];
% 		for j = (length(metrics(i).Ca_mins)-1):-1:1
% 			if isnan(metrics(i).Ca_mins(j))
% 				continue
% 			end
% 			if isnan(metrics(i).Ca_mins(j+1))
% 				continue
% 			end
% 			these_spikes = metrics(i).spike_times(metrics(i).spike_times > metrics(i).Ca_mins(j) & metrics(i).spike_times < metrics(i).Ca_mins(j+1));
% 			this_duration(j) = these_spikes(end) - these_spikes(1);
% 		end
% 		this_duration = 1e-3*mean(nonzeros(this_duration));


% 		if i == 1
% 			% PD neuron

% 			% measure PD neuron duration (1st spike -- last spike)
% 			C =  C + basin_cost(PD_burst_duration_range,this_duration);
% 		end

% 		if i == 2
% 			% LP neuron

% 			% measure LP neuron duration (1st spike -- last spike)
% 			C =  C + basin_cost(LP_burst_duration_range,this_duration);
% 		end

% 		if i == 3
% 			% PY neuron

% 			% measure PY neuron duration (1st spike -- last spike)
% 			C =  C + basin_cost(PY_burst_duration_range,this_duration);
% 		end

% 	end  % end neuron-specific loop

% 	PD_ends = metrics(1).last_spike_loc;
% 	rm_this = PD_ends == 0;
% 	PD_ends = metrics(1).Ca_peaks(:) + PD_ends(:);
% 	PD_ends(rm_this) = [];

% 	PD_starts = metrics(1).first_spike_loc;
% 	rm_this = PD_starts == 0;
% 	PD_starts = metrics(1).Ca_peaks(:) + PD_starts(:);
% 	PD_starts(rm_this) = [];


% 	LP_starts = metrics(2).first_spike_loc;
% 	rm_this = LP_starts == 0;
% 	LP_starts = metrics(2).Ca_peaks(:) + LP_starts(:);
% 	LP_starts(rm_this) = [];

% 	LP_ends = metrics(2).last_spike_loc;
% 	rm_this = LP_ends == 0;
% 	LP_ends = metrics(2).Ca_peaks(:) + LP_ends(:);
% 	LP_ends(rm_this) = [];


% 	PY_starts = metrics(3).first_spike_loc;
% 	rm_this = PY_starts == 0;
% 	PY_starts = metrics(3).Ca_peaks(:) + PY_starts(:);
% 	PY_starts(rm_this) = [];

% 	PY_ends = metrics(3).last_spike_loc;
% 	rm_this = PY_ends == 0;
% 	PY_ends = metrics(3).Ca_peaks(:) + PY_ends(:);
% 	PY_ends(rm_this) = [];

% catch
% 	C = fail_cost;
% 	return
% end


% if nargout == 0
% 	% mark the burst ends and starts
% 	plot(ax(1),[PD_starts PD_starts],[-80 60],'g','LineWidth',3)
% 	plot(ax(1),[PD_ends PD_ends],[-80 60],'r','LineWidth',3)

% 	plot(ax(2),[LP_starts LP_starts],[-80 60],'g','LineWidth',3)
% 	plot(ax(2),[LP_ends LP_ends],[-80 60],'r','LineWidth',3)

% 	plot(ax(3),[PY_starts PY_starts],[-80 60],'g','LineWidth',3)
% 	plot(ax(3),[PY_ends PY_ends],[-80 60],'r','LineWidth',3)
% end

% try

% % measure gaps from PD to LP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% all_gaps = NaN*PD_ends;

% for i = 2:length(all_gaps)-1
% 	[~,idx] = min(abs(PD_ends(i) - LP_starts));
% 	all_gaps(i) = (LP_starts(idx) - PD_ends(i))*1e-3;
% end

% C = C + basin_cost(Gap_PD_LP_range, nanmean(all_gaps));

% % measure gaps from LP to PY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% all_gaps = NaN*LP_ends;

% for i = 2:length(all_gaps)-1
% 	[~,idx] = min(abs(LP_ends(i) - PY_starts));
% 	all_gaps(i) = (PY_starts(idx) - LP_ends(i))*1e-3;
% end

% C = C + basin_cost(Gap_LP_PY_range, nanmean(all_gaps));

% ;;;;;;;;  ;;;;;;;; ;;          ;;;    ;;    ;;  ;;;;;;  
% ;;     ;; ;;       ;;         ;; ;;    ;;  ;;  ;;    ;; 
% ;;     ;; ;;       ;;        ;;   ;;    ;;;;   ;;       
% ;;     ;; ;;;;;;   ;;       ;;     ;;    ;;     ;;;;;;  
% ;;     ;; ;;       ;;       ;;;;;;;;;    ;;          ;; 
% ;;     ;; ;;       ;;       ;;     ;;    ;;    ;;    ;; 
% ;;;;;;;;  ;;;;;;;; ;;;;;;;; ;;     ;;    ;;     ;;;;;;  

% % measure delay from PD -> LP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% all_delays = NaN*PD_starts;

% for i = 2:length(all_delays)-1
% 	[~,idx] = min(abs(PD_starts(i) - LP_starts));
% 	all_delays(i) = (LP_starts(idx) - PD_starts(i))*1e-3;
% end

% C = C + basin_cost(Delay_PD_LP_range, nanmean(all_delays));

% % measure delay from PD -> PY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% all_delays = NaN*PD_starts;

% for i = 2:length(all_delays)-1
% 	[~,idx] = min(abs(PD_starts(i) - PY_starts));
% 	all_delays(i) = (PY_starts(idx) - PD_starts(i))*1e-3;
% end

% C = C + basin_cost(Delay_PD_PY_range, nanmean(all_delays));


% ;;;;;;;;  ;;     ;; ;;;;;;;; ;;    ;; 
% ;;     ;; ;;     ;;    ;;     ;;  ;;  
% ;;     ;; ;;     ;;    ;;      ;;;;   
% ;;     ;; ;;     ;;    ;;       ;;    
% ;;     ;; ;;     ;;    ;;       ;;    
% ;;     ;; ;;     ;;    ;;       ;;    
% ;;;;;;;;   ;;;;;;;     ;;       ;;    

%  ;;;;;;  ;;    ;;  ;;;;;;  ;;       ;;;;;;;; 
% ;;    ;;  ;;  ;;  ;;    ;; ;;       ;;       
% ;;         ;;;;   ;;       ;;       ;;       
% ;;          ;;    ;;       ;;       ;;;;;;   
% ;;          ;;    ;;       ;;       ;;       
% ;;    ;;    ;;    ;;    ;; ;;       ;;       
%  ;;;;;;     ;;     ;;;;;;  ;;;;;;;; ;;;;;;;; 

% C = C + basin_cost(PD_duty_cycle_range,metrics(1).burst_metrics(9));
% C = C + basin_cost(LP_duty_cycle_range,metrics(2).burst_metrics(9));
% C = C + basin_cost(PY_duty_cycle_range,metrics(3).burst_metrics(9));


% % penalize depolarization block in LP

% for i = 1:length(LP_starts)
% 	these_spikes = metrics(2).spike_times(metrics(2).spike_times >= LP_starts(i) & (metrics(2).spike_times <= LP_ends(i)));

% 	allowed_isi_range = [0 3*(these_spikes(end) - these_spikes(1))/length(these_spikes)];

% 	C = C + basin_cost(allowed_isi_range,max(diff(these_spikes)));

% end

% % In PD
% for i = 1:length(PD_starts)
% 	these_spikes = metrics(2).spike_times(metrics(2).spike_times >= PD_starts(i) & (metrics(2).spike_times <= PD_ends(i)));

% 	allowed_isi_range = [0 3*(these_spikes(end) - these_spikes(1))/length(these_spikes)];

% 	C = C + basin_cost(allowed_isi_range,max(diff(these_spikes)));

% end




% catch

% 	C = fail_cost;
% 	return
% end


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


% function c = basin_cost(allowed_range, actual_value)

% 	c = 0;
% 	if actual_value < allowed_range(1)
% 		c = (allowed_range(2) - actual_value)^2/(allowed_range(2) - allowed_range(1));
% 		return
% 	elseif  actual_value > allowed_range(2)
% 		c = (allowed_range(2) - actual_value)^2/(allowed_range(2) - allowed_range(1));
% 		return
% 	else
% 		return 
% 	end


% end % end basin_cost

end % end function 