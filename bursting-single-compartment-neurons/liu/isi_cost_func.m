
function C = isi_cost_func(isis)





C = 1e5;


;;;;;;;;     ;;;    ;;;;;;;;     ;;;    ;;     ;;  ;;;;;;  
;;     ;;   ;; ;;   ;;     ;;   ;; ;;   ;;;   ;;; ;;    ;; 
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;  ;;;; ;;;; ;;       
;;;;;;;;  ;;     ;; ;;;;;;;;  ;;     ;; ;; ;;; ;;  ;;;;;;  
;;        ;;;;;;;;; ;;   ;;   ;;;;;;;;; ;;     ;;       ;; 
;;        ;;     ;; ;;    ;;  ;;     ;; ;;     ;; ;;    ;; 
;;        ;;     ;; ;;     ;; ;;     ;; ;;     ;;  ;;;;;;  


isi_range = [.02 .08];
ibi_range = [.5 1.5 ]; % seconds

n_spikes_per_burst_range = [4 10];


% [V,Ca] = x.integrate;
% cutoff = floor(10e3/x.dt);
% V = V(cutoff:end,1);
% Ca = Ca(cutoff:end,1);



% spike_times = psychopomp.findNSpikes(V,1e3);


%isis = diff(nonnans(spike_times))*1e-3*x.dt;

if isempty(isis)
	C = 1e6;
	return
end


frac_isi = mean(isis>isi_range(1) & isis< isi_range(2));
frac_ibi = mean(isis>ibi_range(1) & isis< ibi_range(2));


% both fractions should sum to 1
C = ((1 - frac_isi + frac_ibi)*100)^2;


% the ratio should be in the range of n_spikes_per_burst
C = C + 1e6*bin_cost(n_spikes_per_burst_range,frac_isi/frac_ibi);

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