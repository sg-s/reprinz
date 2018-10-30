% show the distributions of gbars

x = make_neuron;
x.md5hash;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

channels = x.AB.find('conductance');

for i = 1:7
	subplot(2,4,i); hold on
	hist(all_g(i,:),100)
	title(channels{i})
end


show_these = shuffle(find(~isnan(all_cost)));

for i = 1:length(show_these)


	x.set('*gbar',all_g(:,show_these(i)))
	bursting_cost_func(x)

	pause(2)

	close(gcf)
	drawnow
	close all
end
