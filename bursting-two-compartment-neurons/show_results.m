x = make2C;
show_these = shuffle(find(all_cost == 0));

for i = 1:length(show_these)

	pause(2)

	x.set(x.find('Neurite*gbar'),all_g(:,show_these(i)));	
	two_comp_cost_func(x)

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
	x.set(x.find('Soma*gbar'),all_g(:,show_these(i)));
	x.Soma.NaV.gbar = 0;
	V = x.integrate;
	time = (1:length(V))*x.dt*1e-3;
	plot(ax(i),time,V(:,2),'k')
	set(ax(i),'YLim',[-80 -10],'XLim',[0 3])
end

prettyFig('plw',1);



% for all these neurons, measure the cost without NaV, and the spike height 
cost_wo_NaV = NaN*all_cost;
spike_height = NaN*all_cost;
for i = 1:length(all_cost)
	if rand > .9
		textbar(i,length(all_cost))
	end
	x.set(x.find('Neurite*gbar'),all_g(:,i));	
	x.set(x.find('Soma*gbar'),all_g(:,i));
	x.Soma.NaV.gbar = 0;
	[cost_wo_NaV(i),~,spike_height(i)] = two_comp_cost_func(x);
end

% now t-SNE the gbars and color them by bursting w/o NaV in soma 
R = mctsne(all_g);

figure('outerposition',[200 200 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(R(1,isnan(spike_height)),R(2,isnan(spike_height)),'ko')
plot(R(1,~isnan(spike_height)),R(2,~isnan(spike_height)),'ro')

% now color the t-SNE rep by the spike amplitude 
c = [parula(100); 0 0 0];
colors = zeros(length(all_g),3);
sc = spike_height;
sc = sc - nanmin(sc);
sc = 1+(sc/nanmax(sc))*99;
sc(isnan(sc)) = 101;
colors = c(round(sc),:);

all_g = all_g./all_g(7,:);
R = mctsne(all_g(1:6,:));

figure('outerposition',[200 200 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
scatter(R(1,:),R(2,:),128,colors,'filled')
