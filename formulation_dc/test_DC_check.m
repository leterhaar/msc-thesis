% the feasible solution must be good
assert(DC_check(xstar, dc, wind) == 0);

% an infeasible solution must be not good
assert(DC_check(random_x, dc, wind) == 1);