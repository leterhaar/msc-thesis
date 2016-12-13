% Main function to generate tests
function tests = testScenConsDelta
    tests = functiontests(localfunctions);
end

% Setup function
function setupOnce(testCase) 
    % addpaths if necessary
    if not(exist('DC_f', 'file'))
        addpath('../formulation_dc');
        addpath('../wind');
        addpath('../networks');
    end
    
    % create network and wind models
    testCase.TestData.N = 2;
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

function teardownOnce(testCase)
    
end