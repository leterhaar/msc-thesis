
N_box = 3*ac.N_b;
N_single = 2*ac.N_G;
dim_g = 2*N_box+N_single;

%% test with sdp vars and no j_des

[g, Ws] = AC_g(x, ac, wind.slice(1), 1);

% should return an sdpvar of the right dimensions
assert(isa(g, 'sdpvar'), 'g is not an sdp var');
assert(all(size(g) == [dim_g, 1]), 'g is not the right dimension');
assert(isa(Ws, 'sdpvar'), 'Ws is not a sdpvar');
assert(all(size(Ws) == [2*ac.N_b 2*ac.N_b]), 'Ws not right size');

%% test with sdpvars and a j_des

% pick random box constraint
j_box = randi(N_box);
[g1, Ws] = AC_g(x, ac, wind.slice(1), 1, j_box);

% check the upper limit
g2 = AC_g(x, ac, wind.slice(1), 1, j_box+N_box);

% pick random single constraint
g3 = AC_g(x, ac, wind.slice(1), 1, randi([2*N_box+1 dim_g]));

% check dimensions and types
assert(isa(g1, 'sdpvar'), 'g is not a sdp var');
assert(isa(g2, 'sdpvar'), 'g is not a sdp var');
assert(isa(g3, 'sdpvar'), 'g is not a sdp var');
assert(all(size(g1) == [1, 1]), 'g is not the right dimension');
assert(all(size(g2) == [1, 1]), 'g is not the right dimension');
assert(all(size(g3) == [1, 1]), 'g is not the right dimension');
assert(isa(Ws, 'sdpvar'), 'Ws is not a sdpvar');
assert(all(size(Ws) == [2*ac.N_b 2*ac.N_b]), 'Ws not right size');



%% test with actual values and no j_des
[g, Ws] = AC_g(random_x, ac, wind.slice(1), 1);

% should return a double 
assert(isa(g, 'double'), 'g is not a double');
% without nans
assert(not(any(isnan(g))), 'g has nans');
% of the right size
assert(all(size(g) == [dim_g, 1]), 'g is not the right dimension');

assert(isa(Ws, 'double'), 'Ws is not a sdpvar');
assert(all(size(Ws) == [2*ac.N_b 2*ac.N_b]), 'Ws not right size');


%% test with random values and j_des
% pick a random constraint

for j_box = 1:dim_g
    % check the lower limit
    [g, Ws] = AC_g(random_x, ac, wind.slice(1), 1, j_box);

    % check dimensions and types
    assert(isa(g, 'double'), 'g is not a double');
    assert(all(size(g) == [1, 1]), 'g is not the right dimension');
    assert(isa(Ws, 'double'), 'Ws is not a sdpvar');
    assert(all(size(Ws) == [2*ac.N_b 2*ac.N_b]), 'Ws not right size');
end

