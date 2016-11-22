%% test without j (all constraints)
Ncons = 6*ac.N_b + 2*ac.N_G + 1;
C_all = AC_cons_scen(x, ac, wind.slice(1), 1);

% test length and type
assert(length(C_all) == Ncons, 'constraints not right size');
assert(isa(C_all, 'lmi'), 'constraints not right type');

%% test for every j that a single constraint is returned
C_det = AC_cons_det(x, ac, wind, t);
C_scen = [];
for j = [1:Ncons-1 -1]
    the_constraint = AC_cons_scen(x, ac, wind.slice(1), 1, j);
    assert(length(the_constraint) == 1);
    assert(isa(the_constraint, 'constraint') || isa(the_constraint, 'lmi'));
    C_scen = [C_scen, the_constraint];
end

% test if solution with all constraints individually is still the same as
% all at once
diagnostics2 = optimize([C_det, C_all], Obj, sdpsettings('verbose', 0));
assert(not(diagnostics2.problem), 'Problem optimizing');
xstar2 = values_cell(x);
diagnostics3 = optimize([C_det, C_scen], Obj, sdpsettings('verbose', 0));
assert(not(diagnostics3.problem), 'Problem optimizing');
xstar3 = values_cell(x);

% check if xstar2 is the same as xstar3
assert(all_close(xstar2, xstar3), 'Should be close');