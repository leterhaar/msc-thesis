% Main function to generate tests
function tests = testAACC
    tests = functiontests(localfunctions);
end

% Setup function
function setupOnce(testCase) 
    if not(exist('create_SVM', 'file'))
        addpath('../formulation_SVM');
    end
    
    % create SVM problem
    m = 100;
    d = 10;
    testCase.TestData.svm = create_SVM(d, m);
end

function testLessConstraintsThenACC(testCase)
    ops = sdpsettings('solver', 'mosek');
    svm = testCase.TestData.svm;
    n_agents = 10;
    diam = 4;
    G = random_graph(n_agents, diam, 'rand');
    
    % run the AACC algorithm for only 3 iterations
    [~, agents_AACC] = ACCA(svm.B, svm.delta, svm.deltas, svm.f, ...
                       svm.cons_delta, ...
                       'opt_settings', ops,...
                       'max_its', 6,...
                       'n_agents', n_agents,...
                       'connectivity', G,...
                       'diameter', diam,...
                       'debug', 1);

    % run the ACC algorithm for only 3 iterations
    [~, agents_ACC] = AACC(svm.B, svm.delta, svm.deltas, svm.f, ...
                       svm.cons_delta, ...
                       'opt_settings', ops,...
                       'max_its', 6,...
                       'n_agents', n_agents,...
                       'connectivity', G,...
                       'diameter', diam,...
                       'debug', 1);

    % loop over agents
    for i = 1:n_agents
        for k = 2:6
            % check that the number of deltas is the same or less
            num_ACC = agents_ACC(i).iterations(k).info.num_cons;
            num_AACC = agents_AACC(i).iterations(k).info.num_cons;
            
            verifyLessThan(testCase, num_AACC, num_ACC, sprintf('%i', k));
        end
    end
    
end

function teardownOnce(testCase)
    
end