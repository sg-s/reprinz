% this cost function can be used to fit a two compartment 
% model to data using metrics gleaned from the data


function C = metricsCost(x, data)

x.reset;

x.closed_loop = true;
x.dt = .1;
x.t_end = 10e3;
x.integrate;

V = x.integrate;

C = 0;

% measure metrics in Axon
metrics = xtools.V2metrics(V(:,1),'sampling_rate',10);

C = C + xfit.binCost([data.burst_period - 100, data.burst_period + 100],metrics.burst_period);
C = C + xfit.binCost([data.duty_cycle - .05, data.duty_cycle + .05],metrics.duty_cycle_mean);

% measure minimum and maximum in soma
C = C + 5*xfit.binCost([data.min_V - 2, data.min_V + 2],min(V(:,2)));
C = C + 5*xfit.binCost([data.max_V - 2, data.max_V + 2],max(V(:,2)));