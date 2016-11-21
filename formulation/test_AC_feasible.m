% formulate objective
Obj = AC_f(x, ac, wind, t);

% formulate constraints
C_det = AC_cons_det(x, ac, wind, t);
C_scens = [];
for i = 1:N
    C_scens = [C_scens; AC_cons_scen(x, ac, wind.slice(i), t)];
end

C_all = [C_det, C_scens];

% solve
opt = sdpsettings('verbose', 0, 'solver', 'mosek');
diagnostics = optimize(C_all, Obj, opt);

assert(not(diagnostics.problem), sprintf('Problem solving: %s', diagnostics.info));
assert(value(Obj) > 0, 'Negative objective function does not make sense');

xstar = values_cell(x);

% check feasibility

for i = 1:N
    residuals = AC_g(xstar, ac ,wind.slice(i), t);
    assert(all(residuals < 1e-6), 'infeasible solution');
end