% test the objective function

Obj_test = AC_f(x, ac, wind, 1);

% check dimensions and type
assert(isa(Obj_test, 'sdpvar'), 'Not the right type');
assert(all(size(Obj_test) == [1 1]))

% check for actual values
Obj_test = AC_f(random_x, ac, wind, 1);
assert(isa(Obj_test, 'double'));
assert(all(size(Obj_test) == [1 1]));