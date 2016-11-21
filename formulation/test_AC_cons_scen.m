%% test without j (all constraints)
Ncons = 6*ac.N_b + 2*ac.N_G + 1;
C = AC_cons_scen(x, ac, wind.slice(1), 1);

% test length and type
assert(length(C) == Ncons, 'constraints not right size');
assert(isa(C, 'lmi'), 'constraints not right type');

%% test for every j that a single constraint is returned

for j = [-1 1:Ncons-1]
    C = AC_cons_scen(x, ac, wind.slice(1), 1, j);
    assert(length(C) == 1);
    assert(isa(C, 'constraint') || isa(C, 'lmi'));
end