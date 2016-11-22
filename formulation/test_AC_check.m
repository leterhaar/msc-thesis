%% a feasible solution should not produce an error


% check flag
flag = AC_check(xstar, ac, wind, t);
assert(flag == 0, 'Should be good');

%% a not psd W should give an error
n = 2*ac.N_b;
not_psd_W = ones(n); % start with a PSD matrix

% create a not-psd matrix
while is_psd(not_psd_W)
    not_psd_W = randn(n);
end

% get response
flag = AC_check({not_psd_W, xstar{2:4}}, ac, wind, t);
assert(flag == 1, 'Wf is not PSD');

%% a not rank one Wf should give an error

flag = AC_check({eye(n), xstar{2:4}}, ac, wind, t);
assert(flag == 3, 'Wf is not rank 1');

%% a not rank one Ws should give an error

flag = AC_check({xstar{1}, 100*eye(n), xstar{3:4}}, ac, wind, t);
assert(flag == 4, 'Ws is not rank 1');
flag = AC_check({xstar{1:2}, 100*eye(n), xstar{4}}, ac, wind, t);
assert(flag == 4, 'Ws is not rank 1');

%% a solution that does not satisfy the constraints should give an error

% test with negative values for R (infeasible for deterministic constr)
flag = AC_check({xstar{1}, zeros(n), zeros(n), -ones(2*ac.N_G,1)}, ac, wind, t);
assert(flag == 5, 'Should be infeasible');

% test with zeros for R (infeasible for scenarios)
flag = AC_check({xstar{1:3}, zeros(2*ac.N_G, 1)}, ac, wind, t);
assert(flag == 6, 'Should be infeasible');
