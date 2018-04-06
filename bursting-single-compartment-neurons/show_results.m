x = make2C;
show_these = shuffle(find(all_cost == 0));

for i = 1:length(show_these)

	pause(2)

	x.set(x.find('Neurite*gbar'),all_g(:,show_these(i)));	
	two_comp_cost_func(x,true)

	title(show_these(i))

	drawnow
	

end

% show nice ones
show_these = [175 537 42 592 491 736 934 410 407 914];

figure('outerposition',[0 0 800 900],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:10
	ax(i) = subplot(5,2,i); hold on
	
	ylabel(ax(i),'V_m (mV)')
end

for i = 1:10
	x.set(x.find('Neurite*gbar'),all_g(:,show_these(i)));	
	two_comp_cost_func(x,true,ax(i));
	set(ax(i),'YLim',[-80 -10],'XLim',[0 3])
end

prettyFig('plw',1);

