% create objective function with sdpvar

Obj_test = DC_f(x, dc, wind);

assert(is_scalar(Obj_test), 'wrong dims');
assert(isa(Obj_test, 'sdpvar'), 'wrong type');

% create objective function with doubles

Obj_test = DC_f(random_x, dc, wind);

assert(is_scalar(Obj_test), 'wrong dims');
assert(isa(Obj_test, 'double'), 'wrong type');
assert(not(isnan(Obj_test)), 'has nans');
