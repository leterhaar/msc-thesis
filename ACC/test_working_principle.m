%% test working principle for SVM problem
clear
yalmip('clear');
addpath('../formulation_SVM/');

d = 10;
m = 200;
tol = 1e-2;
svm = create_SVM(d, m);
ops = sdpsettings('verbose', 0);
status = optimize(svm.cons, svm.f(svm.B), ops);
assert(not(status.problem), status.info);
Bstar = value(svm.B);

% retrieve only active constraints
residuals = svm.residuals(Bstar);
active_cons = svm.cons(abs(residuals) < tol);
fprintf('SVM: Removed %d constraints (%2.0f%%)\n', length(residuals)-length(active_cons), (length(residuals)-length(active_cons))/length(residuals)*100);

% optimize again, using only active constraints
status = optimize(active_cons, svm.f(svm.B), ops);
assert(not(status.problem), status.info);
Bstar2 = value(svm.B);

% compare
verify(all_close(Bstar, Bstar2), 'Not close for SVM');

%% Test working principle for a DC problem

yalmip('clear');
N = 5;
dc = DC_model('case14a');
x_sdp = sdpvar(5*dc.N_G,1);
wind = wind_model(dc, 24, 0.2);
wind.dummy(N);
t_wind = 1;
Obj = DC_f_obj(x_sdp, dc, wind, t_wind);
C = DC_f_0(x_sdp, dc, wind, t_wind);

% build scenario constraints
for i = 1:N
    C = [C, DC_f_ineq(x_sdp, i, dc, wind, t_wind)];
end

% solve centralized problem
status = optimize(C, Obj, ops);
assert(not(status.problem), status.info);
xstar = value(x_sdp);

% retrieve only active constraints
assign(x_sdp, xstar);
residuals = check(C);
active_cons = C(abs(residuals) < tol);
fprintf('DC: Removed %d constraints (%2.0f%%)\n', length(residuals)-length(active_cons), (length(residuals)-length(active_cons))/length(residuals)*100);

% optimize again, using only active constraints
status = optimize(active_cons, Obj, ops);
assert(not(status.problem), status.info);
xstar2 = value(x_sdp);

% compare
verify(all_close(xstar2, xstar), 'Not close for DC OPF');

%% Test working principle for AC problem
yalmip('clear');

if not(exist('AC_f', 'file'))
    addpath('../formulation');
    addpath('../wind');
    addpath('../misc');
    addpath('../networks');
end

% load models
N = 50;      % number of scenarios used for optimization
t = 1; % timestep used for this demonstration (todo: add for everything)

% load network and wind models
ac = AC_model('case14a');
ac.set_WPG_bus(9);
wind = wind_model(ac, 24, 0.2);

% generate a number of scenarios
wind.generate(N);

% optimization settings
ops = sdpsettings('solver', 'mosek', 'verbose', 0);

% create SDPvars for all variables
Wf = sdpvar(2*ac.N_b);
Wu = sdpvar(2*ac.N_b);
Wd = sdpvar(2*ac.N_b);
R = sdpvar(2*ac.N_G, 1);
x_cell = {Wf, Wu, Wd, R};

% define objective function
Obj = AC_f(x_cell, ac, wind, t);

% define general constraints
C = AC_cons_det(x_cell, ac, wind, t);

for i = 1:N
    C = [C, AC_cons_scen(x_cell, ac, wind.slice(i), t)];
end
status = optimize(C, Obj, ops);
assert(not(status.problem), status.info);
xstar = values_cell(x_cell);

% retrieve only active constraints
assign_cell(x_cell, xstar);
residuals = check(C);
active_cons = C(abs(residuals) < tol);
fprintf('AC: Removed %d constraints (%2.0f%%)\n', length(residuals)-length(active_cons), (length(residuals)-length(active_cons))/length(residuals)*100);

% optimize again, using only active contraints
status = optimize(active_cons, Obj, ops);
assert(not(status.problem), status.info);
xstar2 = values_cell(x_cell);

% compare
[close, maxdiff] = all_close(xstar, xstar2, 1e-3);
verify(close, 'Not close: %g', maxdiff);