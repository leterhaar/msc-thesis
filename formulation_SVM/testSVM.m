%% Main function to generate tests
function tests = testSVM
    tests = functiontests(localfunctions);
end

%% Test Functions
function setupOnce(testCase)
    d = 50; m = 200;
    svm = create_SVM(d, m);
    one_class = svm.ys > 0;

    initfig('test SVM', 3);
    plot(svm.xs(1, one_class), svm.xs(2, one_class), 'x');
    plot(svm.xs(1, ~one_class), svm.xs(2, ~one_class), 'o');
    theaxis = axis;
    drawline = [-100 100];
    plot(svm.Bstar(2)*drawline, -svm.Bstar(1)*drawline, '--');
    axis(theaxis);
    testCase.TestData.svm = svm;
    testCase.TestData.d = d;
    testCase.TestData.m = m;
end

function testOptimalSolutionsClose(testCase)
    svm = testCase.TestData.svm;
    optimize(svm.cons, svm.f(svm.B), sdpsettings('verbose', 0));
    verifyEqual(testCase, svm.Bstar, value(svm.B), 'RelTol', 1e-10);

end

function testResiduals(testCase)
    svm = testCase.TestData.svm;
    d = testCase.TestData.d;
    m = testCase.TestData.m;
    newB = svm.Bstar;
    newB(randperm(d,ceil(d/10))) = rand(ceil(d/10), 1);

    assign(svm.B, newB);
    residuals1 = check(svm.cons);

    residuals2 = nan(m,1);
    for i = 1:m
        residuals2(i) = svm.residual(newB, i);
    end

    residuals3 = svm.residuals(newB);

    verifyEqual(testCase, residuals1, residuals2, 'RelTol', 1e-10);
    verifyEqual(testCase, residuals1, residuals3', 'RelTol', 1e-10);
end

function teardownOnce(testCase)  % do not change function name
    close(3);
end