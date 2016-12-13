% create SVM problem
addpath('../formulation_SVM');

opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');

% run algorithm
[xstar, agents] = ACC(svm.B, svm.delta, svm.deltas, svm.f, svm.cons_delta, ...
                   'verbose', 1,...
                   'debug', 1, ...
                   'opt_settings', opt_settings,...
                   'n_agents', m, ...
                   'diameter', 2);
%% calculate convergence and feasibility
m = length(agents);
K = length(agents(1).iterations);

convergence = nan(K,m);
feasibility = nan(K,m);
time_per_iteration = nan(K,m);
optimal_objective = svm.f(svm.Bstar);
% p = progress('Checking constraints', m);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence(k,i) = abs(agents(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(svm.B, agents(i).iterations(k).x)
        feasibility(k,i) = sum(check(svm.cons) < -1e-6) / N_svm * 100;
        
        % store times
        time_per_iteration(k,i) = agents(i).iterations(k).time;
    end
%     p.ping();
end

%% plot
initfig('ACC iterations', 1);
ax = subplot(211, 'YScale', 'log');
grid on
hold on
plot(convergence);
ylabel('|f(x_k^i) - f(x^*) |')
title(sprintf('ACC convergence for SVM with d=%i, m=%i', d, N_svm));

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility);
ylabel('% violated');
xlabel('iterations');

initfig('ACC timing', 2);
plot(time_per_iteration);

assert(all_close(xstar, svm.Bstar, 1e-4), 'Not close');
assert(sum(feasibility(k,:)) == 0, 'Not all feasible');

%% run the same but then with residual function instead of check
% (must be faster)


% run algorithm
[xstar2, agents2] = ACC(svm.B, svm.delta, svm.deltas, svm.f, svm.cons_delta, ...
                   'verbose', 1,...
                   'debug', 1, ...
                   'opt_settings', opt_settings,...
                   'n_agents', m, ...
                   'diameter', 2, ...
                   'residuals', svm.residual_delta);

%% calculate convergence and feasibility
m = length(agents2);
K = length(agents2(1).iterations);

convergence = nan(K,m);
feasibility = nan(K,m);
time_per_iteration = nan(K,m);
optimal_objective = svm.f(svm.Bstar);
% p = progress('Checking constraints', m);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence(k,i) = abs(agents2(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(svm.B, agents2(i).iterations(k).x)
        feasibility(k,i) = sum(check(svm.cons) < -1e-6) / N_svm * 100;
        
        % store times
        time_per_iteration(k,i) = agents2(i).iterations(k).time;
    end
%     p.ping();
end

%% plot
initfig('ACC iterations', 3);
ax = subplot(211, 'YScale', 'log');
grid on
hold on
plot(convergence);
ylabel('|f(x_k^i) - f(x^*) |')
title(sprintf('ACC convergence for SVM with d=%i, m=%i fn', d, N_svm));

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility);
ylabel('% violated');
xlabel('iterations');

initfig('ACC timing', 4);
plot(time_per_iteration);

assert(all_close(xstar2, svm.Bstar, 1e-4), 'Not close');
assert(sum(feasibility(k,:)) == 0, 'Not all feasible');
