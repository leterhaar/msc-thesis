% test equivalence between optimizer and optimize command for svm problem

% create SVM problem
m = 100;
d = 20;
svm = create_SVM(d,m);
opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');

Obj = svm.f(svm.B);

% create constraints and test equivalence with 5 random Bs
for i = 1:5
    assign(svm.B, rand(d,1));
    
    res1 = check(svm.cons);
    res2 = [];
    for j = 1:m
        assign(svm.delta, svm.deltas(j, :))
        res2 = [res2; check(svm.cons_delta)];
    end
    assert(all_close(res1, res2));
end

% solve optimization problem using optimize
status = optimize(svm.cons, Obj, opt_settings);
assert(not(status.problem), status.info);
assert(all(check(svm.cons) > -1e-6), 'Not feasible')
xstar1 = value(svm.B);

% solve optimization problem using optimizer
solvert = optimizer(svm.cons_delta, Obj, opt_settings, svm.delta, svm.B);
merged = [];
for i = 1:m
    merged = [merged, solvert(svm.deltas(i, :), 'nosolve')];
end
[xstar2, problem, msg] = merged([]);
assert(not(problem), msg);
assign(svm.B, xstar2);
for i = 1:m
    assign(svm.delta, svm.deltas(i, :))
    assert(all(check(svm.cons_delta) > -1e-6), 'Not feasible');
end

assert(all_close(xstar1, xstar2), 'Not the same');
assert(all_close(xstar1, svm.Bstar), 'Not optimal');