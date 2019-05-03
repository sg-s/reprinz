
function C = exact_fit_cost_func(x,data)


step_1_cost = 1e4;
step_2_cost = 1e3;
step_3_cost = 1e2;
C = step_1_cost + step_2_cost;

 ;;;;;;  ;;;;;;;; ;;;;;;;; ;;;;;;;;        ;;   
;;    ;;    ;;    ;;       ;;     ;;     ;;;;   
;;          ;;    ;;       ;;     ;;       ;;   
 ;;;;;;     ;;    ;;;;;;   ;;;;;;;;        ;;   
      ;;    ;;    ;;       ;;              ;;   
;;    ;;    ;;    ;;       ;;              ;;   
 ;;;;;;     ;;    ;;;;;;;; ;;            ;;;;;; 


% step 1 -- fit slow wave using metrics
tolerance = 2;
burst_peak_range = [data.mu_peak - tolerance*data.sigma_peak; data.mu_peak + tolerance*data.sigma_peak];
burst_min_range = [data.mu_min - tolerance*data.sigma_min; data.mu_min + tolerance*data.sigma_min];
burst_period_range = [data.mu_period - tolerance*data.sigma_period; data.mu_period + tolerance*data.sigma_period];

% copy conductances densities from cell body to neurite 
x.Neurite.ACurrent.gbar = x.CellBody.ACurrent.gbar;
x.Neurite.CaS.gbar = x.CellBody.CaS.gbar;
x.Neurite.CaT.gbar = x.CellBody.CaT.gbar;
x.Neurite.HCurrent.gbar = x.CellBody.HCurrent.gbar;
x.Neurite.KCa.gbar = x.CellBody.KCa.gbar;
x.Neurite.Kd.gbar = x.CellBody.Kd.gbar;
x.Neurite.Leak.gbar = x.CellBody.Leak.gbar;


% turn off 2nd compartment
x.synapses(1).resistivity = Inf;
x.synapses(2).resistivity = Inf;

x.sim_dt = .1;
x.dt = .1;



x.t_end = 20e3;
V = x.integrate;
Vf(:,1) = sgolayfilt(V(:,1),1,1001);
V = V(10e4:end,1);
Vf =  Vf(10e4:end,1);

V = Vf;

% find peaks 
[ons,offs] = computeOnsOffs(V > median(V));

if length(ons) < 5
	return
end

burst_peaks = NaN*ons;
burst_peak_locs = NaN*ons;
burst_trough_locs = NaN*ons;
burst_mins = NaN*ons;

for i = 2:length(ons)-1
	[burst_peaks(i), temp] = max(V(ons(i-1):ons(i)));
	burst_peak_locs(i) = temp + ons(i-1);
	[burst_mins(i), temp] = min(V(ons(i-1):ons(i)));
	burst_trough_locs(i) = temp + ons(i-1);
end

burst_peak_locs = burst_peak_locs(2:end-1);
burst_periods = diff(burst_peak_locs);
burst_peaks = burst_peaks(2:end-1);
burst_mins = burst_mins(2:end-1);

% cost for peaks
burst_peak_cost = bin_cost(burst_peak_range,mean(burst_peaks))*1e4/3;
burst_min_cost = bin_cost(burst_min_range,mean(burst_mins))*1e4/3;
burst_period_cost = bin_cost(burst_period_range,mean(burst_periods))*1e4/3;

this_step_cost = burst_peak_cost + burst_min_cost + burst_period_cost;

if this_step_cost > 0
	C = this_step_cost + step_2_cost;
	return
end

 ;;;;;;  ;;;;;;;; ;;;;;;;; ;;;;;;;;      ;;;;;;;  
;;    ;;    ;;    ;;       ;;     ;;    ;;     ;; 
;;          ;;    ;;       ;;     ;;           ;; 
 ;;;;;;     ;;    ;;;;;;   ;;;;;;;;      ;;;;;;;  
      ;;    ;;    ;;       ;;           ;;        
;;    ;;    ;;    ;;       ;;           ;;        
 ;;;;;;     ;;    ;;;;;;;; ;;           ;;;;;;;;; 


% now we try to match it exactly to the slow wave

a = find(abs(V - mean(burst_mins)) < .1,1,'first');
V = V(a+1:a+9e4);


this_step_cost = step_2_cost*(1- rsquare(V,data.slow_wave));
C = this_step_cost;

if this_step_cost > step_3_cost

	return
end



 ;;;;;;  ;;;;;;;; ;;;;;;;; ;;;;;;;;      ;;;;;;;  
;;    ;;    ;;    ;;       ;;     ;;    ;;     ;; 
;;          ;;    ;;       ;;     ;;           ;; 
 ;;;;;;     ;;    ;;;;;;   ;;;;;;;;      ;;;;;;;  
      ;;    ;;    ;;       ;;                  ;; 
;;    ;;    ;;    ;;       ;;           ;;     ;; 
 ;;;;;;     ;;    ;;;;;;;; ;;            ;;;;;;;  






if nargout == 0
	figure, hold on
	plot(data.slow_wave,'k')
	plot(V,'r')
end




	function c = bin_cost(allowed_range,actual_value)

		if isnan(actual_value)
			actual_value = allowed_range(2)*10;
		end

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