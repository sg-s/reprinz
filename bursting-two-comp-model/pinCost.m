% pin a neuron to a clamping voltage
% and measure the cost of that pin

function C = pinCost(x,data)



x.reset;

x.closed_loop = false;
x.dt = .1;
x.sim_dt = .1;
x.t_end = 9e3;

x.CellBody.V = data.V0(1);

x.V_clamp = [NaN*data.V0 data.V0];


V = x.integrate;


C = sum((V(:,2)).^2);


x.V_clamp = NaN*x.V_clamp;
x.closed_loop = true;
x.integrate;
V = x.integrate;

% also compare voltage traces directly
C = C + 2*xtools.voltageCost(data.V0,V(:,2),200);


if isnan(C)
	C = 1e3;
	return
end