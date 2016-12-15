clear
yalmip('clear');

% addpaths if necessary
if not(exist('DC_f', 'file'))
    addpath('../formulation_dc');
    addpath('../wind');
    addpath('../networks');
end
    
% create network and wind models
N = 10;
dc = DC_model('case14a');
dc.set_WPG_bus(9);
wind = wind_model(dc, 24, 0.2);
wind.generate(N);

opt_settings = sdpsettings('solver', 'gurobi', 'verbose', 0);
n_agents = 3;
diam = 2;
x_sdp_full = sdpvar(5*dc.N_G, 24);
delta_sdp = sdpvar(1, 3*24);
f = @(x) DC_f(x, dc, wind);
constraints = DC_cons_scen_delta(x_sdp_full, dc, delta_sdp);
deltas = [wind.P_w', ...
          max(0, wind.P_m'),...
          max(0, -wind.P_m')];
default_constraint = DC_cons_det(x_sdp_full, dc, wind);
residuals = @(x, delta, j) DC_g_delta(x, dc, delta, j);

%% run ACC algorithm

[xstar, agents] = ACC(x_sdp_full, delta_sdp, deltas, f, constraints, ...
                      'verbose', 1,...
                      'opt_settings', opt_settings,...
                      'default_constraint', default_constraint,...
                      'diameter', diam,...
                      'n_agents', n_agents,...
                      'debug', 1,...
                      'residuals', residuals,...
                      'use_selector', true);

%% run centralized problem
C_scens = [];
for i = 1:N
    C_scens = [C_scens, DC_cons_scen_single(x_sdp_full, dc, wind.slice(i))];
end

status = optimize([default_constraint, C_scens], f(x_sdp_full), opt_settings);
assert(not(status.problem), status.info);
xstar_centralized = value(x_sdp_full);                  
%% calculate convergence and feasibility
K = length(agents(1).iterations);
C_all = [default_constraint, C_scens];
convergence = nan(K,n_agents);
feasibility = nan(K,n_agents);
time_per_iteration = nan(K,n_agents);
optimal_objective = f(xstar_centralized);
optimizations_run = zeros(K,1);
ACC_active_deltas = nan(K,n_agents);
for i = 1:n_agents
    for k = 1:K
        % calculate difference with centralized objective
        convergence(k,i) = abs(agents(i).iterations(k).J ...
                                                - optimal_objective);

        % calculate feasibility percentage
        assign(x_sdp_full, agents(i).iterations(k).x)
        feasibility(k,i) = sum(check(C_all) < -1e-5) / N * 100;

        % store times
        time_per_iteration(k,i) = agents(i).iterations(k).time;

        % store no of constraints
        ACC_active_deltas(k,i) = ...
                            size(agents(i).iterations(k).active_deltas, 1);
        
        % store total number of iterations run
        if k > 1
            optimizations_run(k) = optimizations_run(k) + ...
                                    agents(i).iterations(k).info.optimized;
        end
    end
end

%% plot
initfig('ACC iterations', 1);
ax = subplot(211, 'YScale', 'log');
semilogy(convergence);
grid on
ylabel('|f(x_k^i) - f(x^*) |')
title('Convergence');

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility);
ylabel('% violated');
xlabel('iterations');

initfig('ACC timing', 2);
plot(time_per_iteration);

initfig('ACC Active Deltas', 3);
plot(ACC_active_deltas);
plot(repmat(length(C_scens), K, 1), '--')



%% Make assertions
assert(all_close(xstar, xstar_centralized), 'Not the same');
assert(sum(feasibility(end, :)) == 0, 'Not all feasible in the end');
