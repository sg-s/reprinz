make_stg
show_these = shuffle(find(~isnan(all_cost)));

for i = 1:length(show_these)

	pause(3)
	close all

	x.set(x.find('*gbar'),all_g(:,show_these(i)));	
	STG_cost_function(x)

	drawnow
	close(gcf)

	

end