C = AC_cons_det(x, ac, wind, 1);

% test length and type
assert(length(C) == 3*ac.N_b+7, 'Det constraints not right size');
assert(isa(C, 'lmi'), 'Det constraints not right type');
