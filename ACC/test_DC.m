dc = DC_model('case14a');
dc.set_WPG_bus(9);
wind_dc = wind_model(dc, 24, 0.2);

x_dc = sdpvar(5*dc.N_G,1);
t_wind = 1;

% formulate problem
Obj = DC_f_obj(x_dc, dc, wind, t_wind);
C_all = DC_f_0(x_dc, dc, wind, t_wind);
for i = 1:N
    C_all = [C_all, DC_f_ineq(x_dc, i, dc, wind, t_wind)];
end
opt = sdpsettings('verbose', 0);
diagnostics = optimize(C_all, Obj, opt);
xstar = value(x_dc);

assert(not(diagnostics.problem), 'Problem optimizing');

% extract active params
active_params = [];
for i = 1:N
    active_params = [active_params; DC_f_check(xstar, i, dc, wind, t_wind)];
end

% build new constraint set with only active params
C_act = DC_f_0(x_dc, dc, wind, t_wind);
for p = active_params'
    i = p(1);
    j = p(2);
    C_act = [C_act, DC_f_ineq(x_dc, i, dc, wind, t_wind, j)];
end

diagnostics2 = optimize(C_act, Obj, opt);
assert(not(diagnostics2.problem), 'Problem optimizing');
xstar2 = value(x_dc);

% check if they are all close
assert(all_close(xstar, xstar2), 'Not close');