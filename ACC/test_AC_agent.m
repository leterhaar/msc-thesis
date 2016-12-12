% initialize agent
ag = AC_agent(ac, ac_wind, 1, 1, 3);

Ncons = 6*ac.N_b + 2*ac.N_G + 1;
% check size of constraints
assert(length(ag.C_1_params) == 3*Ncons, 'not right size');

% check if objective function has value
assert(isa(ag.Obj, 'sdpvar'), 'obj should be sdpvar');
assert(not(isnan(value(ag.Obj))), 'optimization should return value');

% check that active params are stored right
all_params = [];
for i = 1:3
    x = values_cell(ag.x);
    act_params = AC_active(x, ac, ac_wind.slice(i), 1);
    all_params = [all_params; ones(length(act_params),1)*i act_params];
end
assert(all(size(ag.A{1}) == size(all_params)), 'Params not same size');
assert(all(all(ag.A{1} == all_params)), 'Active params not matching');