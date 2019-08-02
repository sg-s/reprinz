% this script shows the best solutions we found
function show_results()

p = twoComp.setup;
data = twoComp.getData;


d = fileparts(mfilename('fullpath'));

% load data results
load([d filesep 'results.xfit'],'-mat');


% remove duplicates
[all_params, idx] = unique(all_params,'rows');
all_cost = all_cost(idx);
all_hash = all_hash(idx);

% check that the cost is valid
current_cost_hash = hashlib.md5hash(which(func2str(p.SimFcn)));

for i = 1:length(all_hash)
	corelib.textbar(i,length(all_hash))
	if ~strcmp(char(all_hash(i)),current_cost_hash)
		% recompute cost
		C = Inf(length(data),1);
		for j = 1:length(data)
			p.data = data(j);
			C(j) = p.evaluate(all_params(i,:));
		end

		all_cost(i) = min(C);
		all_hash(i) = current_cost_hash;

	end
end

save([d filesep 'results.xfit'],'all_hash','all_cost','all_params')




% show some passing ones and throw in the real ones for sport

show_these = [714 579 967 244 1173 760 1161 115 765 445 775 1177 635 945 51 886 837 407 802 413 924 193 1072];

figure('outerposition',[300 300 801 1100],'PaperUnits','points','PaperSize',[801 1100]); hold on

real_data = veclib.shuffle(1:18);
real_data = real_data(1:4)
idx = 1;

% add some noise
p.x.CellBody.add('CurrentNoise','noise_amplitude',.5);

for i = 1:18
	subplot(6,3,i); hold on

	if ismember(i,real_data)
		time = (1:length(data(idx).V0))*1e-4;
		plot(time,data(idx).V0,'k')
		set(gca,'XLim',[0 5])
		idx = idx + 1;
	else

		p.x.set(p.parameter_names,all_params(show_these(i),:))
		p.x.closed_loop = true;
		p.x.reset;
		p.x.integrate;
		p.x.integrate;
		V = p.x.integrate;
		time = (1:length(V))*p.x.dt*1e-3;
		plot(time,V(:,2)+ randn(length(V),1)*5e-2,'k')

		%title([mat2str(show_these(i)) '  cost = ' strlib.oval(all_cost(show_these(i)))])
	end
	axis off
end

figlib.pretty('PlotLineWidth',1)

return















[~,show_these] = sort(all_cost);

figure('outerposition',[300 300 1801 901],'PaperUnits','points','PaperSize',[1200 901]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
	ph(i) = plot(NaN,NaN,'k');
end


for i = 1:length(show_these)
	


	if any(isnan(all_params(show_these(i),:)))
		continue
	end

	idx = rem(i,6) + 1;

	p.x.set(p.parameter_names,all_params(show_these(i),:))
	p.x.closed_loop = true;
	p.x.reset;
	p.x.integrate;
	p.x.integrate;
	V = p.x.integrate;
	time = (1:length(V))*p.x.dt*1e-3;
	ph(idx).XData = time;
	ph(idx).YData = V(:,2) + randn(length(V),1)*5e-2;

	title(ax(idx),[mat2str(show_these(i)) '  cost = ' strlib.oval(all_cost(show_these(i)))])

	pause(1)

end