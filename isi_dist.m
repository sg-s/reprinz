% plot distributions of ISIs

figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on

for i = 1:length(all_spikes)
	isi = diff(all_spikes{i})*1e-4;
	[hy,hx] = histcounts(isi,linspace(0,1,1e3));
	hx = hx(1:end-1) + mean(diff(hx));
	hy = hy/sum(hy);
	plot(hx,hy)

end