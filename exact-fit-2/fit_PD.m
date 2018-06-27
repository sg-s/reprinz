
load('../exact-fit/isolated_PD')
V0 = V;
T = 1e3;
V = filtfilt(ones(T,1),T,V);
V = V(1:9e4);
V0 = V0(1:9e4);
dV = [NaN; diff(V)];
dV0 = [NaN; diff(V0)];

% find burst period
[V_max,loc]=findpeaks(V);
clear data
data.slow_wave = V;
data.V0 = V0;
data.mu_period = mean(diff(loc));
data.sigma_period = std(diff(loc));
data.mu_peak = mean(V_max);
data.sigma_peak = std(V_max);
V_min=findpeaks(-V);
data.mu_min = mean(-V_min);
data.sigma_min = std(-V_min);

% first just fit the slow wave
x = make2C;


p = procrustes('particleswarm');
p.x = x;

p.parameter_names = [x.find('CellBody*gbar'); 'CellBody.Ca'] ;

p.data = data;

M = length(p.parameter_names);



% neuron conductances
%      A   CaS  CaT   H    KCa   Kd    Leak Ca    
ub = [2e3  2e2  125   1    2e3   2e3   10   3e3 ];
lb = [0    0    0     0    0     0     0    .05 ];



p.lb = lb;
p.ub = ub;

p.sim_func = @exact_fit_cost_func;

load('best_seed.mat','S')
p.seed = S;




return





i = 1;

p.seed = [all_g(:,i); .05];



return




figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on

subplot(1,2,1); hold on
plot(V0,'k')
plot(V,'r')
set(gca,'XLim',[1e3 2e4])
subplot(1,2,2); hold on
plot(V,dV)


N = xolotl.findNSpikes(V,-30);
spiketimes = xolotl.findNSpikeTimes(V,N,-30);



x = make2C;


p = procrustes('particleswarm');
p.x = x;

p.data.LeMassonMatrix = procrustes.V2matrix(V,[-80 50],[-20 30]);
p.data.V = V;
p.data.spiketimes = spiketimes;

p.parameter_names = [x.find('Neurite*gbar'); x.find('CellBody*gbar'); 'Neurite.tau_Ca'; 'Neurite.Ca'; 'synapses(1).resistivity'; 'temperature'];

M = length(p.parameter_names);

%         A    CaS  CaT   H    KCa    Kd   Leak      NaV 
g_ub =  [500  60   100    .1    100   2e3   1        0];
g_lb =  [10   10   10     .001  10    100   1e-3     0];

% neuron conductances

p.ub = [g_ub g_ub(1:end-1) 200 3e3 .01   30];
p.lb = [g_lb g_lb(1:end-1) 20  .05 .0001 -10] ;


p.seed = rand(M,1).*p.ub(:); % random seed

p.sim_func = @exact_fit_cost_func;
p.options.MaxTime = 300;




