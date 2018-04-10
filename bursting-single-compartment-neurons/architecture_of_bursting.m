
pHeader;


%% The Architecture of Bursting Neurons
% In this document, I look at the various ways in which we can get a bursting neuron. 

% load the data
load('reprinz_1c_dalek.bio.brandeis.edu.mat')

% trim the data
all_g = all_g(:,find(~isnan(all_cost)));
ratio_g = all_g./all_g(end,:);


%%
% First, I show some sample solutions picked at random. Note that there appears to be some evidence of degenerate solutions. 


x = make_neuron;
x.t_end = 3e3;

channels = x.AB.find('conductance');

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1234)); 
figure('outerposition',[0 0 1000 900],'PaperUnits','points','PaperSize',[1000 900]); hold on
clear ax

for i = 5:-1:1
	ax(i*2-1) = subplot(5,2,i*2-1);
	idx = randi(length(all_g));
	x.set(x.find('*gbar'),all_g(:,idx));

	x.reset;
	x.integrate;
	V = x.integrate;

	time = x.dt*(1:length(V))*1e-3;
	plot(ax(i*2-1),time,V,'k','LineWidth',1)
	ax(i*2-1).Box = 'off';

	% show the conductance profile
	ax(i*2) = subplot(5,2,i*2); hold on
	for j = 1:size(all_g,1)
		stem(j,all_g(j,idx))
	end
	set(ax(i*2),'YScale','log','YLim',[1e-2 1e4],'YTick',[1e-2 1 1e2 1e4])
end

for i = 1:8
	set(ax(i),'XColor','w')
end

ax(10).XTick = 1:7;
ax(10).XTickLabel = channels;
ax(10).XTickLabelRotation = 45;

xlabel(ax(9),'Time (s)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% What does the distribution of individual conductances look like? 

figure('outerposition',[0 0 1000 900],'PaperUnits','points','PaperSize',[1000 900]); hold on

c = lines(7);

for i = 1:7
	subplot(3,3,i); hold on
	[hy,hx] = histcounts(all_g(i,:),linspace(0,max(all_g(i,:))*1.2,100));
	hx = hx(1:end-1) + mean(diff(hx));
	stairs(hx,hy,'Color',c(i,:))
end

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% What do the pairwise distributions look like? Here we plot the ratio of every conductance vs. the sodium conductance. Overall, there appears to be little correlation, suggesting that we are finding a wide range of solutions. 

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on

c = lines(7);

for i = 1:6
	subplot(2,3,i); hold on
	scatter(all_g(i,:),all_g(end,:),32,[.4 .4 .4],'filled','Marker','o','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.9);
	set(gca,'YColor',c(end,:))
	set(gca,'XColor',c(i,:))
	xlabel([ 'g_{' channels{i} '} (uS/mm^2)'])
	ylabel([ 'g_{' channels{end} '} (uS/mm^2)'])
end

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% How can we visualize all these conductances, to see if there is some structure here? I can embed all conductances using t-SNE to see if something sticks out. 

R = mctsne(all_g);

figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
scatter(R(1,:),R(2,:),32,[.4 .4 .4],'filled','Marker','o','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.9);
axis off

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Embedding of ratios
% 


R = mctsne(ratio_g(1:6,:));

figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
scatter(R(1,:),R(2,:),32,[.4 .4 .4],'filled','Marker','o','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.9);
axis off

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


