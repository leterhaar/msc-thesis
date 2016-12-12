%% test to check whether the problem formulation and solution
% with 1 scenario is the same (i.e. whether DC_f_ineq and 
% DC_f_ineq_delta with the same delta yield the same solutions)

% define common problem
x_sdp = sdpvar(dc.N_G*5,1,'full');
delta_sdp = sdpvar(1,2, 'full');
Obj = DC_f_obj(x_sdp, dc, wind.slice(1), t_wind);
Cons_det = DC_f_0(x_sdp, dc, wind.slice(1), t_wind);
opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');


% check whether residuals are the same for 10 random instances of x
Cons_scen_delta = DC_f_ineq_delta(x_sdp, delta_sdp, dc, t_wind);
for i = 1:N
    x = rand(dc.N_G*5, 1);
    test_delta = [wind.P_w(t_wind, i) wind.P_m(t_wind, i)];
    assign(x_sdp, x);
    assign(delta_sdp, test_delta);
    Cons_fixed_delta = DC_f_ineq_delta(x_sdp, test_delta, dc, t_wind);
    Cons_scen_different_i = DC_f_ineq(x_sdp, i, dc, wind, t_wind);
    
    residuals1 = check(Cons_scen_different_i);
    residuals2 = check(Cons_scen_delta);
    residuals3 = check(Cons_fixed_delta);
    
    assert(all_close(residuals1, residuals3), 'Residuals not the same');
    assert(all_close(residuals1, residuals2), 'Residuals not the same')
end

% solve problem with scenario constraint
Cons_scen = DC_f_ineq(x_sdp, 1, dc, wind, t_wind);
status = optimize([Cons_det; Cons_scen], Obj, opt_settings);
assert(not(status.problem), status.info);
xstar1 = value(x_sdp);
value(Obj)

% check feasibility
assert(all(check(Cons_det) > -1e-6), 'Infeasible solution');
assert(all(check(Cons_scen) > -1e-6), 'Infeasible solution');


% solve problem with scenario constraint using optimizer
delta = [wind.P_w(t_wind, 1) wind.P_m(t_wind, 1)];
Cons_scen_delta = DC_f_ineq_delta(x_sdp, delta_sdp, dc, t_wind);

the_solver = optimizer([Cons_det; Cons_scen_delta], Obj, ...
                       opt_settings, delta_sdp, x_sdp);
                   
the_second_solver = optimizer([Cons_det; Cons_scen], Obj, ...
                        opt_settings, [], x_sdp);
[xstar_r, problem, msg] = the_solver(delta);
assert(not(problem), msg);

[xstar_r2, problem, msg] = the_second_solver([]);
assert(not(problem), msg);

assert(all_close(xstar_r, xstar_r2, 1e-3), 'Not close');
assert(all_close(xstar_r, xstar, 1e-3), 'Not close');

% check feasibility
assign(x_sdp, xstar_r);
assign(delta_sdp, delta);
assert(all(check(Cons_det) > -1e-5), 'Deterministic infeasible solution');
assert(all(check(Cons_scen_delta) > -1e-6), 'Scenario infeasible solution');

% solve the problem using the same function but with optimize (w/o -r)
Cons_scen_fixed_delta = DC_f_ineq_delta(x_sdp, delta, dc, t_wind);
status = optimize([Cons_det, Cons_scen_fixed_delta], Obj, opt_settings);
assert(not(status.problem), status.info);
xstar_fixed_delta = value(x_sdp);

% check feasibility
assert(all(check(Cons_det) > -1e-6), 'Infeasible solution');
assert(all(check(Cons_scen_fixed_delta) > -1e-6), 'Infeasible solution');



% check equivalence
assert(all_close(xstar1, xstar_fixed_delta), 'Solutions with functions not the same');
assert(all_close(xstar1, xstar_r, 1e-2), 'Solutions with optimizer not the same');

