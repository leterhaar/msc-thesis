% create SVM problem
clear
yalmip('clear');
addpath('../formulation_SVM');

m_svm = 5;      % number of agents

opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');

% run algorithm
[xstar, agents] = ACC(svm.B, svm.delta, svm.deltas, svm.f, svm.cons_delta, ...
                   'verbose', 1,...
                   'debug', 1, ...
                   'opt_settings', opt_settings,...
                   'n_agents', m_svm);
%% calculate convergence and feasibility
m_svm = length(agents);
K = length(agents(1).iterations);

convergence = nan(K,m_svm);
feasibility = nan(K,m_svm);
time_per_iteration = nan(K,m_svm);
optimal_objective = svm.f(svm.Bstar);
p = progress('Checking constraints', m_svm);
for i = 1:m_svm
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
    p.ping();
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
