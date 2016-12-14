% Test the results for the ACC algorithm for simple DC OPF 
% against centralized solution

% build problem
m = 4; % number of agents
f = @(x) DC_f_obj(x, dc, wind, t_wind);
default_constraint = DC_f_0(x_sdp, dc, wind, t_wind);
constraints_delta = DC_f_ineq_delta(x_sdp, delta_sdp, dc, t_wind);
opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');
%%
% build deltas and scenario constraints
deltas = [];
scenario_constraints = [];
for i = 1:N
    deltas = [deltas; wind.P_w(t_wind, i), ...
                  max(0, -wind.P_m(t_wind, i)), ...
                  max(0, wind.P_m(t_wind, i))];
              
    scenario_constraints = [scenario_constraints, ... 
                            DC_f_ineq(x_sdp, i, dc, wind, t_wind)];
end

C_all = [default_constraint, scenario_constraints];
% solve centralized problem
status = optimize(C_all, f(x_sdp), opt_settings);
assert(not(status.problem), status.info);
xstar_centralized = value(x_sdp);

% solve ACC algorithm
[xstar_acc, agents] = ACC(x_sdp, delta_sdp, deltas, f, constraints_delta, ...
                          'verbose', 1, ...
                          'default_constraint', default_constraint, ...
                          'n_agents', m,...
                          'debug', 1,...
                          'opt_settings', opt_settings,...
                          'max_its', 20);

%% calculate convergence and feasibility
m = length(agents);
K = length(agents(1).iterations);

convergence = nan(K,m);
feasibility = nan(K,m);
time_per_iteration = nan(K,m);
optimal_objective = f(xstar_centralized);
ACC_active_deltas = nan(K,m);
optimizations_run = zeros(K,1);

for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence(k,i) = abs(agents(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(x_sdp, agents(i).iterations(k).x)
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
ax = subplot(211);
hold off
plot(convergence);
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
yyaxis left
plot(time_per_iteration);
ylabel('Time per iteration');
yyaxis right
plot(optimizations_run);
ylabel('Optimizations run');

initfig('ACC Active Deltas', 3);
plot(ACC_active_deltas);
plot(repmat(length(scenario_constraints), K, 1), '--')


%% Make assertions
assert(all_close(xstar_acc, xstar_centralized), 'Not the same');
assert(sum(feasibility(end, :)) == 0, 'Not all feasible in the end');