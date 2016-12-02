%% test without j (all constraints)
Ncons = (4*dc.N_G + 2*dc.N_l)*N_t;
C_all = DC_cons_scen2(x, dc, wind.slice(1));

% test length and type
assert(length(C_all) == Ncons, 'constraints not right size');
assert(isa(C_all, 'lmi'), 'constraints not right type');

%% test for every j that a single constraint is returned
if not(exist('random_j', 'var'))
    random_j = nan(1,10);
    for i = 1:10
        random_j(i) = randi(Ncons);
    end
end

for j = random_j
    the_constraint = DC_cons_scen2(x, dc, wind.slice(1), j);
    assert(length(the_constraint) == 1);
    assert(isa(the_constraint, 'constraint') || isa(the_constraint, 'lmi'));
end

%% test equivalence with faster notation

C2 = DC_cons_scen(x, dc, wind);

info = optimize([C_det, C2], Obj, sdpsettings('verbose', 0));
assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));

assert(all_close(value(x), xstar), 'Not the same');