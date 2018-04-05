% first, assemble all the data

assemble_all_data

% filter data based on cost
keep_these = (find(C == 0));
C = C(keep_these);
M = M(:,keep_these);
g = g(:,keep_these);


burst_duration_range(:,1) = [0.4490  0.7150]; % PD, sec
burst_duration_range(:,2) = [0.2860 .5120]; % LP, sec
burst_duration_range(:,3) = [0.38 .68]; % PY, sec

n_spikes_per_burst_range(:,1) = [4 30];
n_spikes_per_burst_range(:,2) = [4 30];
n_spikes_per_burst_range(:,3) = [4 30];

cycle_period_range = [1.23 1.788]; % mean +/- 1 SD

duty_cycle_range(:,1) = [.3450 .4250];
duty_cycle_range(:,2) = [.2050 .3230];
duty_cycle_range(:,3) = [.2940 .4020];


% define gaps A_end -> B_start
Gap_range(:,1) = [.112 .33];
Gap_range(:,2) = [.112 .33];
Gap_range(:,3) = [NaN NaN];

Delay_range(:,1) = [NaN NaN]; % seconds
Delay_range(:,2) = [.6340 .9720]; % seconds
Delay_range(:,3) = [.9250 1.3570]; % seconds



figure('outerposition',[30 30 2000 1200],'PaperUnits','points','PaperSize',[1200 600]); hold on
clear ax
for i = 1:18
	ax(i) = subplot(3,6,i); hold on
end


% show durations 
for i = 1:3
	this_ax = (i-1)*6 + 1;
	[hy,hx] = histcounts(M(i,:),linspace(0,1,100));
	hx = hx(1:end-1) + mean(diff(hx));
	stairs(ax(this_ax), hx,hy,'k')
	xlabel(ax(this_ax),'Burst duration (s)')

	plot(ax(this_ax),[burst_duration_range(1,i) burst_duration_range(1,i)],[0 max(hy)],'r')
	plot(ax(this_ax),[burst_duration_range(2,i) burst_duration_range(2,i)],[0 max(hy)],'r')
	set(ax(this_ax),'XLim',[0 1])

end

% show # of spikes/burst
for i = 1:3
	this_ax = (i-1)*6 + 2;
	[hy,hx] = histcounts(M(i+3,:),linspace(0,50,50));
	hx = hx(1:end-1) + mean(diff(hx));
	stairs(ax(this_ax), hx,hy,'k')
	xlabel(ax(this_ax),'#spikes/burst')

	plot(ax(this_ax),[n_spikes_per_burst_range(1,i) n_spikes_per_burst_range(1,i)],[0 max(hy)],'r')
	plot(ax(this_ax),[n_spikes_per_burst_range(2,i) n_spikes_per_burst_range(2,i)],[0 max(hy)],'r')
	set(ax(this_ax),'XLim',[0 50])
end


% show burst period 
for i = 1:3
	this_ax = (i-1)*6 + 3;
	[hy,hx] = histcounts(M(i+6,:),linspace(0,3,100));
	hx = hx(1:end-1) + mean(diff(hx));
	stairs(ax(this_ax), hx,hy,'k')
	xlabel(ax(this_ax),'Burst period (s)')

	plot(ax(this_ax),[cycle_period_range(1) cycle_period_range(1)],[0 max(hy)],'r')
	plot(ax(this_ax),[cycle_period_range(2) cycle_period_range(2)],[0 max(hy)],'r')
	set(ax(this_ax),'XLim',[0 3])
end

% show duty cycle
for i = 1:3
	this_ax = (i-1)*6 + 4;
	[hy,hx] = histcounts(M(i+9,:),linspace(0.1,.5,100));
	hx = hx(1:end-1) + mean(diff(hx));
	stairs(ax(this_ax), hx,hy,'k')
	xlabel(ax(this_ax),'Duty cycle)')

	plot(ax(this_ax),[duty_cycle_range(1,i) duty_cycle_range(1,i)],[0 max(hy)],'r')
	plot(ax(this_ax),[duty_cycle_range(2,i) duty_cycle_range(2,i)],[0 max(hy)],'r')
	set(ax(this_ax),'XLim',[0 1])
end


% gaps
for i = 1:3
	this_ax = (i-1)*6 + 5;
	[hy,hx] = histcounts(M(i+12,:),linspace(-.5,.5,100));
	hx = hx(1:end-1) + mean(diff(hx));
	stairs(ax(this_ax), hx,hy,'k')
	xlabel(ax(this_ax),'Gap (s)')

	plot(ax(this_ax),[Gap_range(1,i) Gap_range(1,i)],[0 max(hy)],'r')
	plot(ax(this_ax),[Gap_range(2,i) Gap_range(2,i)],[0 max(hy)],'r')
	set(ax(this_ax),'XLim',[-.5 .5])
end

% delays 
for i = 2:3
	this_ax = (i-1)*6 + 6;
	[hy,hx] = histcounts(M(i+15,:),linspace(.5,1.5,100));
	hx = hx(1:end-1) + mean(diff(hx));
	stairs(ax(this_ax), hx,hy,'k')
	xlabel(ax(this_ax),'Delay (s)')

	plot(ax(this_ax),[Delay_range(1,i) Delay_range(1,i)],[0 max(hy)],'r')
	plot(ax(this_ax),[Delay_range(2,i) Delay_range(2,i)],[0 max(hy)],'r')
	set(ax(this_ax),'XLim',[.5 1.5])

end

prettyFig;