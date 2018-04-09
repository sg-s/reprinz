% show the distributions of gbars

x = make_neuron;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

channels = x.AB.find('conductance');

for i = 1:7
	subplot(2,4,i); hold on
	hist(all_g(i,:),100)
	title(channels{i})
end
