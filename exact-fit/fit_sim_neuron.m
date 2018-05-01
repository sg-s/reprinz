



% conversion from Prinz to phi
A = 0.0628;
vol = A; % this can be anything, doesn't matter
f = 14.96; % uM/nA
tau_Ca = 200;
phi = (2*f*96485*vol)/tau_Ca;

channels = {'NaV','CaT','CaS','ACurrent','KCa','Kd','HCurrent'};
prefix = 'prinz-approx/';
gbar = [1000 25  60 500  50  1000 .1];
E =    [50   30  30 -80 -80 -80   -20];

x = xolotl;
x.add('AB','compartment','Cm',10,'A',A,'vol',vol,'phi',phi,'Ca_out',3000,'Ca_in',0.05,'tau_Ca',tau_Ca);




for i = 1:length(channels)
	x.AB.add([prefix channels{i}],'gbar',gbar(i),'E',E(i));
end

g0 = x.get(x.find('*gbar'));

x.t_end = 5e3;

x.integrate; 
V = x.integrate;

p = procrustes('particleswarm');
p.x = x;
p.V_clamp = [V];

p.parameter_names = [x.find('AB*gbar'); 'AB.Ca'];

M = length(p.parameter_names);


% neuron conductances
p.ub = [1e3  300  300  10    100   3e3   3e3 x.AB.Ca_out];
p.lb = [.1   .1   .1   1e-3   .1   500   1e2  x.AB.Ca_in];


p.seed = rand(M,1).*p.ub(:); % random seed

x.set(p.parameter_names,p.seed);

I_clamp_before = x.integrate([],p.V_clamp');

p.sim_func = @exact_fit_cost_func;
p.options.MaxTime = 300;
p.fit;

x.set(p.parameter_names,p.seed)
V = x.integrate;

figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
subplot(2,3,1:2); hold on
time = (1:length(V))*x.dt*(1e-3);
plot(time,p.V_clamp,'k')
title('Behavior of "data"')


subplot(2,3,4:5); hold on
time = (1:length(V))*x.dt*(1e-3);
plot(time,V,'r')
title('Behavior of model fit')

g = x.get(x.find('*gbar'));

subplot(2,3,3); hold on
plot(g0,g,'o')
set(gca,'XScale','log','YScale','log','XTick',[1e-2 1 1e2 1e4],'YTick',[1e-2 1 1e2 1e4])
xlabel('Original conductances')
xlabel('Fit conductances')
prettyFig();



figure('outerposition',[300 300 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
subplot(2,1,1); hold on
plot(time,I_clamp_before,'k');

subplot(2,1,2); hold on
I_clamp = x.integrate([],p.V_clamp');
plot(time,I_clamp,'r')