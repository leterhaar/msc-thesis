%% init test vars for DC

dc = DC_model('case14a');
dc.set_WPG_bus(9);
N = 10;
wind = wind_model(dc, 24, 0.2);
wind.dummy(10);
N_t = 24;

x = sdpvar(5*dc.N_G,N_t);
random_x = rand(5*dc.N_G,N_t);

Obj = DC_f(x, dc, wind);
C = DC_cons_det(x, dc, wind);

for i = 1:N
    C = [C, DC_cons_scen(x, dc, wind.slice(i))];
end

diagnostics = optimize(C, Obj, sdpsettings('verbose', 0));
assert(not(diagnostics.problem), 'Problem optimizing');
xstar = value(x);

test_sequence = {...
    'DC_f',...
    'DC_cons_det',...
    'DC_g',...
    'DC_cons_scen',...
    };