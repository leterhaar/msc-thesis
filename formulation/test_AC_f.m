% test the objective function

Obj = AC_f(x, ac, wind, 1);

% check dimensions and type
assert(isa(Obj, 'sdpvar'), 'Not the right type');
assert(all(size(Obj) == [1 1]))

% check for actual values
Obj = AC_f(random_x, ac, wind, 1);
assert(isa(Obj, 'double'));
assert(all(size(Obj) == [1 1]));