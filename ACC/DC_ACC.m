clear
yalmip('clear');
% add helper path
addpath('../experiment');
N = 10;
m = 3;
% init experiment
init_experiment(...
    'model_name', 'case14a', ... adapted 14 bus network
    'model_formulation', 'DC', ...
    'wind_N', N, ... scenarios
    'wind_Nm', ceil(N/m), ... scenarios per agent
    'algorithm', 'ACC', ...
    'wind_dummy', 1);
t_wind = 8;
% generate random connection graph with fixed diameter
dm = 2;

% build problem
x_sdp = sdpvar(5*dc.N_G,1, 'full');
delta_sdp = sdpvar(1,2, 'full');
f = @(x) DC_f_obj(x, dc, wind, t_wind);
default_constraint = DC_f_0(x_sdp, dc, wind, t_wind);
constraints_delta = DC_f_ineq_delta(x_sdp, delta_sdp, dc, t_wind);
opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');
%%
% build deltas
deltas = [wind.P_w(t_wind, :)' wind.P_m(t_wind, :)'];
[xstar_acc, agents] = ACC(x_sdp, delta_sdp, deltas, f, constraints_delta, ...
                          'verbose', 1, ...
                          'default_constraint', default_constraint, ...
                          'n_agents', m,...
                          'diameter', dm,...
                          'debug', 1,...
                          'opt_settings', opt_settings,...
                          'max_its', 10);