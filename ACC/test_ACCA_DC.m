
diam = 2;
m = 4; % number of agents
G = random_graph(m, diam, 'rand');
f = @(x) DC_f_obj(x, dc, wind, t_wind);
default_constraint = DC_f_0(x_sdp, dc, wind, t_wind);
constraints_delta = @(x, delta) DC_f_ineq_delta(x, delta, dc, t_wind);
ops = sdpsettings('verbose', 0, 'solver', 'gurobi');
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
lmis_delta = DC_f_ineq_delta(x_sdp, delta_sdp, dc, t_wind);


C_all = [default_constraint, scenario_constraints];
% solve centralized problem
status = optimize(C_all, f(x_sdp), ops);
assert(not(status.problem), status.info);
xstar_centralized = value(x_sdp);

% run the ACCA algorithm
[xstar_ACCA, agents_ACCA] = ACCA_fcn2(x_sdp, deltas, f, ...
                            constraints_delta, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'max_its', 20,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 1,...
                           'default_constraint', default_constraint);

% run the ACC algorithm
[xstar_ACC, agents_ACC] = ACC(x_sdp, delta_sdp, deltas, f, ...
                            lmis_delta, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'max_its', 20,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 1,...
                           'default_constraint', default_constraint);
                       
%% calculate convergence and feasibility
K = length(agents_ACCA(1).iterations);
convergence_ACCA = nan(K,m);
feasibility_ACCA = nan(K,m);
time_per_iteration_ACCA = nan(K,m);
optimal_objective = f(xstar_centralized);
optimizations_run_ACCA = zeros(K,1);
no_cons_used_ACCA = nan(K,m);
p = progress('Checking constraints', m);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACCA(k,i) = abs(agents_ACCA(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(x_sdp, agents_ACCA(i).iterations(k).x)
        feasibility_ACCA(k,i) = sum(check(C_all) < -1e-6) / length(C_all) * 100;
        
        % store times
        time_per_iteration_ACCA(k,i) = agents_ACCA(i).iterations(k).time;
        
        % store total number of iterations run
        if k > 1
            optimizations_run_ACCA(k) = optimizations_run_ACCA(k) + ...
                                    agents_ACCA(i).iterations(k).info.optimized;
            
            no_cons_used_ACCA(k,i) = agents_ACCA(i).iterations(k).info.num_cons;

        end
    end
    p.ping();
end

%% calculate convergence and feasibility
K = length(agents_ACC(1).iterations);

convergence_ACC = nan(K,m);
feasibility_ACC = nan(K,m);
time_per_iteration_ACC = nan(K,m);
optimizations_run_ACC = zeros(K,1);
no_cons_used_ACC = nan(K,m);
p = progress('Checking constraints', m);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACC(k,i) = abs(agents_ACC(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(x_sdp, agents_ACC(i).iterations(k).x)
        feasibility_ACC(k,i) = sum(check(C_all) < -1e-6) / length(C_all) * 100;
        
        % store times
        time_per_iteration_ACC(k,i) = agents_ACC(i).iterations(k).time;
        
        % store total number of iterations run
        if k > 1
            optimizations_run_ACC(k) = optimizations_run_ACC(k) + ...
                                agents_ACC(i).iterations(k).info.optimized;
            no_cons_used_ACC(k,i) = agents_ACC(i).iterations(k).info.num_cons;
        end
    end
    p.ping();
end

%% plot
initfig('ACC iterations', 1);
ax = subplot(211, 'YScale', 'log');
grid on
hold on
plot(convergence_ACCA, 'color', green);
plot(convergence_ACC, ':', 'color', blue);
ylabel('|f(x_k^i) - f(x^*) |')
title(sprintf('ACC convergence for DC OPF with N=%i ACCA', N));

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility_ACCA, 'color', green);
plot(feasibility_ACC, ':', 'color', blue);
ylabel('% violated');
xlabel('iterations');

initfig('ACC timing', 2);
ax = subplot(211);
grid on
hold on
plot(time_per_iteration_ACCA, 'color', green);
plot(time_per_iteration_ACC, ':', 'color', blue);
ylabel('Time per iteration');

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(optimizations_run_ACCA, 'color', green);
plot(optimizations_run_ACC, ':', 'color', blue);
ylabel('Optimizations run');
xlabel('Iteration');

initfig('No of constraints used', 3);
plot(no_cons_used_ACCA, 'color', green);
plot(no_cons_used_ACC, ':', 'color', blue);

figure(1);
%% make assertions
verify(all_close(xstar_ACC, svm.Bstar, 1e-4), 'ACC not optimal');
verify(all_close(xstar_ACCA, svm.Bstar, 1e-4), 'ACCA not optimal');
verify(sum(feasibility_ACC(k,:)) == 0, 'ACC not all feasible');
verify(sum(feasibility_ACCA(k,:)) == 0, 'ACCA not all feasible');