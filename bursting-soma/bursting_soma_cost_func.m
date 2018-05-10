% cost function that will be used to tune 
% a mulit-compartment neuron model


function C = bursting_soma_cost_func(x,~,~)

;;;;;;;;     ;;;    ;;;;;;;;     ;;;    ;;     ;;  ;;;;;;  
;;     ;;   ;; ;;   ;;     ;;   ;; ;;   ;;;   ;;; ;;    ;; 
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;  ;;;; ;;;; ;;       
;;;;;;;;  ;;     ;; ;;;;;;;;  ;;     ;; ;; ;;; ;;  ;;;;;;  
;;        ;;;;;;;;; ;;   ;;   ;;;;;;;;; ;;     ;;       ;; 
;;        ;;     ;; ;;    ;;  ;;     ;; ;;     ;; ;;    ;; 
;;        ;;     ;; ;;     ;; ;;     ;; ;;     ;;  ;;;;;;  


slow_wave_range = [20 30];
cycle_period_range = [.8 2]; 
duty_cycle_range = [.3 .5];

% integrate, etc. 
[V,~,~,I] = x.integrate;
x.t_end = 15e3;
cutoff = floor(5e3/x.dt);
V = V(cutoff:end,1);
I = I(cutoff:end,:);
T = length(V)*x.dt*1e-3; % seconds
time = (1:length(V))*x.dt*1e-3; % seconds

% find peaks 
[ons,offs] = computeOnsOffs(V > median(V));

if length(ons) < 5
	C = 1e3;
	return
end

burst_peaks = NaN*ons;
burst_peak_locs = NaN*ons;
burst_trough_locs = NaN*ons;
burst_troughs = NaN*ons;

for i = 2:length(ons)-1
	[burst_peaks(i), temp] = max(V(ons(i-1):ons(i)));
	burst_peak_locs(i) = temp + ons(i-1);
	[burst_troughs(i), temp] = min(V(ons(i-1):ons(i)));
	burst_trough_locs(i) = temp + ons(i-1);
end

burst_peak_locs = burst_peak_locs(2:end-1);
burst_peaks = burst_peaks(2:end-1);
burst_troughs = burst_troughs(2:end-1);

% difure out up and down swings
dV = [0; diff(V)];

% Hcurrent should contribute to upswing 
curr_contrib = xolotl.contributingCurrents(V,I);
H_contrib = sum((curr_contrib == 4).*(dV>0))/sum(dV>0);


% KCa should contribute to downswing 
KCa_contrib = sum((curr_contrib == 5).*(dV<0))/sum(dV<0);

contrib_cost = 2*bin_cost([.1 .9],H_contrib) + bin_cost([.1 .9],KCa_contrib);


period_cost = 2*bin_cost(cycle_period_range,mean(diff(burst_peak_locs*x.dt*1e-3)));



amplitude_cost = bin_cost(slow_wave_range,mean(burst_peaks - burst_troughs));

% duty cycle
V_midpoint = (mean(burst_peaks) + mean(burst_troughs))/2;
duty_cycle = mean(V>V_midpoint);
dc_cost = bin_cost(duty_cycle_range,duty_cycle);



C = dc_cost + amplitude_cost + period_cost + contrib_cost;

if nargout == 0
	figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
	plot(time,V,'k')

	disp(['Period cost = ' oval(period_cost)])
	disp(['Duty cycle cost = ' oval(dc_cost)])
	disp(['H current conribution is = ' oval(H_contrib)])
	disp(['KCa current conribution is = ' oval(KCa_contrib)])
	disp(['Amplitude  cost = ' oval(amplitude_cost)])
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