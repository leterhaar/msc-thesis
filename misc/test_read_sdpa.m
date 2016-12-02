% A non-existent file should return an error
try 
    read_sdpa('foobar');
    doesnotexist = 0;
catch
    doesnotexist = 1;
end
assert(logical(doesnotexist), 'non existing file should throw error');

% run with existing file, check if matrices are correct
problem = read_sdpa('sdpa_test_file.dat');
C = problem.C;
As = problem.As;
bs = problem.bs;
X = problem.X;


real_bs = [10 20];
assert(all_close(bs, real_bs));

real_C = -[1 0 0 0; 0 2 0 0; 0 0 3 0; 0 0 0 4];
assert(all_close(C, real_C));

real_As = {[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]; ...
      [0 0 0 0; 0 1 0 0; 0 0 5 2; 0 0 2 6]};
assert(all_close(As, real_As));

assert(isa(X, 'sdpvar'));
assert(all(size(X) == [4 4]));

%% test it with a real problem

problem = read_sdpa('../sdplib/hinf1.dat-s');
C = problem.C;
As = problem.As;
bs = problem.bs;
X = problem.X; 

Obj = trace(C * X);
Cons = [X >= 0];
for i = 1:length(As)
    Cons = [Cons, trace(As{i}*X) == bs(i)];
end

diagnostics = optimize(Cons, Obj, sdpsettings('verbose', 0, 'solver', 'mosek'));
assert(not(diagnostics.problem), 'Failed: %s', diagnostics.info)
assert(abs(value(Obj)+2.0326) < 1e-3)