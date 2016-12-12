d = 100; m = 2000;
svm = create_SVM(d, m);

one_class = svm.ys > 0;

initfig('test SVM');
plot(svm.xs(1, one_class), svm.xs(2, one_class), 'x');
plot(svm.xs(1, ~one_class), svm.xs(2, ~one_class), 'o');
theaxis = axis;
drawline = [-100 100];
plot(svm.Bstar(2)*drawline, -svm.Bstar(1)*drawline, '--');
axis(theaxis);

optimize(svm.cons, svm.f(svm.B), sdpsettings('verbose', 0));
value(svm.B);

assert(all_close(svm.Bstar, value(svm.B)));

% make a partial violating B
newB = svm.Bstar;
newB(randperm(d,ceil(d/10))) = rand(ceil(d/10), 1);
tic
assign(svm.B, newB);
residuals1 = check(svm.cons);
toc
tic
residuals2 = nan(m,1);
for i = 1:m
    residuals2(i) = svm.residual(newB, i);
end
toc
tic 
residuals3 = svm.residuals(newB);
toc
assert(all_close(residuals1, residuals2));
assert(all_close(residuals1, residuals3'));
