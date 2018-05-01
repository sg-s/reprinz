x = make2C;
x.sha1hash;
show_these = shuffle(find(all_cost == 0));

for i = 1:length(show_these)

	if spike_height(show_these(i)) < 5
		continue
	end

	if isnan(spike_height(show_these(i)))
		continue
	end

	pause(1)
	x.set(x.find('Neurite*gbar'),all_g(:,show_these(i)));	
	x.set(x.find('Soma*gbar'),all_g(:,show_these(i)));
	%x.Soma.NaV.gbar = 0;
	x.integrate;
	V = x.integrate;
	time = (1:length(V))*x.dt*1e-3;
	figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	plot(time,V(:,2),'k')
	set(gca,'YLim',[-80 -10],'XLim',[0 3])
	title(show_these(i))
	drawnow
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
