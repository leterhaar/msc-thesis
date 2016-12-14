% Main function to generate tests
function tests = testScenConsDelta
    tests = functiontests(localfunctions);
end

% Setup function
function setupOnce(testCase) 
    yalmip('clear');
    % addpaths if necessary
    if not(exist('DC_model', 'file'))
        addpath('../wind');
        addpath('../networks');
    end
    
    % create network and wind models
    testCase.TestData.N = 10;
    dc = DC_model('case14a');
    dc.set_WPG_bus(9);
    wind = wind_model(dc, 24, 0.2);
    wind.generate(testCase.TestData.N);
    
    % store in testcase
    testCase.TestData.dc = dc;
    testCase.TestData.wind = wind;
end

function testSameResiduals(testCase)
    % restore data
    dc = testCase.TestData.dc;
    wind = testCase.TestData.wind;
    N = testCase.TestData.N;
    randX = rand(5*dc.N_G, 24);
    x_sdp = sdpvar(5*dc.N_G, 24);
    assign(x_sdp, randX);
    
    % test for every scenario that the residuals are the same
    for i = 1:N
        % formulate with wind model and check
        cons_scen = DC_cons_scen_single(x_sdp, dc, wind.slice(i));
        residuals1 = check(cons_scen);
        
        % formulate deltas and check
        delta = [   wind.P_w(:,i)' ...
                    max(0, wind.P_m(:,i))' ...
                    max(0, -wind.P_m(:,i))'];
        cons_scen_delta = DC_cons_scen_delta(x_sdp, dc, delta);
        residuals2 = check(cons_scen_delta);
        
        verifyEqual(testCase, residuals1, residuals2, 'RelTol', 1e-10);
        
        % calculate residuals directly
        residuals3 = DC_g_delta(randX, dc, delta);
        verifyEqual(testCase, residuals1, residuals3, 'RelTol', 1e-10);
    end
end

function testUsingSelector(testCase)
    % restore data
    dc = testCase.TestData.dc;
    wind = testCase.TestData.wind;
    randX = rand(5*dc.N_G, 24);
    
    % formulate deltas and check
    i = 1;
        delta = [   wind.P_w(:,i)' ...
                    max(0, wind.P_m(:,i))' ...
                    max(0, -wind.P_m(:,i))'];
    all_residuals = DC_g_delta(randX, dc, delta);
    Ncons = length(all_residuals);
    all_residuals2 = nan(Ncons, 1);
    for j = 1:Ncons
        all_residuals2(j) = DC_g_delta(randX, dc, delta, j);
    end
    
    verifyEqual(testCase, all_residuals2, all_residuals);
    
end

function testSolvableWithSliceOfScenarios(testCase)
    
    % test if using just 1 or 2 scenarios, the problem is solvable
    dc = testCase.TestData.dc;
    wind = testCase.TestData.wind;
    N = testCase.TestData.N;
    x_sdp = sdpvar(5*dc.N_G, 24);
    C_det = DC_cons_det(x_sdp, dc, wind);
    Obj = DC_f(x_sdp, dc, wind);
    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    for m = 1:N
        cons = C_det;
        for i = 1:m
            delta = [   wind.P_w(:,i)' ...
                    max(0, wind.P_m(:,i))' ...
                    max(0, -wind.P_m(:,i))'];
            cons = [cons, DC_cons_scen_delta(x_sdp, dc, delta)];
        end
        status = optimize(cons, Obj, ops);
        verifyTrue(testCase, not(status.problem), ...
                    sprintf(['With %i scens and delta: ' status.info], m));
        xstar1 = value(x_sdp);
        
        cons = C_det;
        for i = 1:m
            cons = [cons, DC_cons_scen_single(x_sdp, dc, wind.slice(i))];
        end
        status = optimize(cons, Obj, ops);
        verifyTrue(testCase, not(status.problem), ...
                    sprintf(['With %i scens and wind: ' status.info], m));
        xstar2 = value(x_sdp);
        verifyEqual(testCase, xstar1, xstar2, 'RelTol', 1e-5);
    end
end

function testSameSolution(testCase)

    % restore data
    dc = testCase.TestData.dc;
    wind = testCase.TestData.wind;
    N = testCase.TestData.N;
    x_sdp = sdpvar(5*dc.N_G, 24);
    
    C_det = DC_cons_det(x_sdp, dc, wind);
    Obj = DC_f(x_sdp, dc, wind);
    ops = sdpsettings('verbose', 0, 'solver', 'mosek');
    
    C_scens = [];
    deltas = nan(N, 72);
    for i = 1:N
        C_scens = [C_scens, DC_cons_scen(x_sdp, dc, wind.slice(i))];
        deltas(i,:) = [   wind.P_w(:,i)' ...
                    max(0, wind.P_m(:,i))' ...
                    max(0, -wind.P_m(:,i))'];
    end
    
    deltas2 = [wind.P_w', ...
          max(0, wind.P_m'),...
          max(0, -wind.P_m')];
      
    verifyEqual(testCase, deltas, deltas2);
    
    status = optimize([C_det, C_scens], Obj, ops);
    verifyTrue(testCase, not(status.problem), status.info);
    xstar1 = value(x_sdp);
    
    % check with delta notation
    C_scens2 = [];
    for i = 1:N
        C_scens2 = [C_scens2, DC_cons_scen_delta(x_sdp, dc, deltas(i, :))];
    end
    status = optimize([C_det, C_scens2], Obj, ops);
    verifyTrue(testCase, not(status.problem), status.info);
    xstar2 = value(x_sdp);
    
    verifyEqual(testCase, xstar1, xstar2, 'RelTol', 1e-5);
    
    % test the same with optimizer
    delta_sdp = sdpvar(1,72);
    C_scen3 = DC_cons_scen_delta(x_sdp, dc, delta_sdp);
    solvert = optimizer([C_det, C_scen3], Obj, ops, delta_sdp, x_sdp);
    merged = [];
    for i = 1:N
        merged = [merged; solvert(deltas(i, :), 'nosolve')];
    end
    
    [xstar3, problem, msg] = merged();
    verifyTrue(testCase, not(problem), msg{:});
    
    verifyEqual(testCase, xstar3, xstar1, 'RelTol', 1e-4);
    
end

function testOldvsNew(testCase)

    

end



function teardownOnce(testCase)
    
end