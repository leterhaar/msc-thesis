Ncons = N_t*7;

C_det = DC_cons_det(x, dc, wind);

assert(length(C_det) == Ncons, 'wrong size');
assert(isa(C_det, 'lmi'), 'wrong type');