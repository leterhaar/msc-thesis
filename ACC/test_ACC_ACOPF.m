% initialize problem

clc
clear
yalmip('clear');

if not(exist('AC_f', 'file'))
    addpath('../formulation');
    addpath('../wind');
    addpath('../misc');
    addpath('../networks');
end

%% load models
N_t = 24;           % optimization horizon
N = 20;             % number of scenarios used for optimization
t = 1;              % timestep used for this demonstration (todo: add for everything)
tol = 1e-4;

% load network and wind models
ac = AC_model('case14a');
ac.set_WPG_bus(9);
wind = wind_model(ac, N_t, 0.2);

% generate a number of scenarios
wind.generate(N);
% wind.use_extremes(t);
% N = 2;

% optimization settings
ops = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug', 1);

% connectivity matrix
n_agents = 5;
G = ones(n_agents)-diag(ones(n_agents,1));
diam = 1;
%% create objective and constraint functions

% create SDPvars for all variables
Wf = sdpvar(2*ac.N_b);
Wu = sdpvar(2*ac.N_b);
Wd = sdpvar(2*ac.N_b);
R = sdpvar(2*ac.N_G, 1);
x_cell = {Wf, Wu, Wd, R};

% define objective function
objective_fcn = @(x) AC_f(x, ac, wind, t);

% define general constraints
default_constraints = AC_cons_det(x_cell, ac, wind, t);

% define deltas
deltas = [wind.P_w(t,:)' wind.P_m(t,:)'];

% define constraint function
cons_fcn = @(x, delta, j_des) AC_cons_scen(x, ac, ... fake wind object using delta
                    struct('P_w', delta(1), 'P_m', delta(2)), 1, j_des);

% define residual function
residuals = @(x, delta, j_des) -AC_g(x, ac, ... again, fake wind object for delta
                    struct('P_w', delta(1), 'P_m', delta(2)), 1, j_des);


%% solve centralized problem
C_all = default_constraints;
C_all_delta = default_constraints;
for i = 1:N
    C_all = [C_all, AC_cons_scen(x_cell, ac, wind.slice(i), t)];
    
    C_all_delta = [C_all_delta, cons_fcn(x_cell, deltas(i, :), 0)];
end
Ncons = (length(C_all) - length(default_constraints)) / N;
%%
tic
status = optimize(C_all, objective_fcn(x_cell), ops);
toc
verify(not(status.problem), status.info);
xstar_cent = values_cell(x_cell);
% 
% status = optimize(C_all_delta, objective_fcn(x_cell), ops);
% verify(not(status.problem), status.info);
% xstar_cent_using_delta = values_cell(x_cell);
% 
% verify(all_close(xstar_cent, xstar_cent_using_delta), ...
%                                             'Cons fcn is not equivalent');
%% run algorithm

[xstar, agents_ACC] = ACC_fcn_cell(x_cell, deltas, objective_fcn, cons_fcn, ...
                             'default_constraint', default_constraints,...
                             'residuals', residuals,...
                             'use_selector', true,...
                             'opt_settings', ops,...
                             'connectivity', G,...
                             'tolerance', tol,...
                             'diameter', diam,...
                             'max_its', 10,...
                             'n_agents', n_agents,...
                             'verbose', 1,...
                             'debug', 1,...
                             'x0', []);
%%                         
verify(AC_check(xstar, ac, wind, t) == 0, 'Infeasible solution');

%% calculate convergence and feasibility
K = length(agents_ACC(1).iterations);
convergence_ACC = nan(K,n_agents);
feasibility_ACC = nan(K,n_agents);
time_per_iteration_ACC = nan(K,n_agents);
optimal_objective = objective_fcn(xstar_cent);
optimizations_run_ACC = zeros(K,1);
no_cons_used_ACC = nan(K,n_agents);
no_cons_active = nan(K,n_agents);
p = progress('Checking constraints', n_agents * K);
timing_ACC = zeros(K,1);

for i = 1:n_agents
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACC(k,i) = abs(agents_ACC(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        infeasible = 0;
        for s = 1:N
            residuals = -AC_g(agents_ACC(i).iterations(k).x, ac, wind.slice(s), t);
            infeasible = infeasible + (sum(residuals < -tol) / length(residuals) / N * 100);
        end
        feasibility_ACC(k,i) = infeasible;
        
        % store times
        time_per_iteration_ACC(k,i) = agents_ACC(i).iterations(k).time;
        timing_ACC(k) = timing_ACC(k) + agents_ACC(i).iterations(k).time;

        
        % store total number of iterations run
        if true
            optimizations_run_ACC(k) = optimizations_run_ACC(k) + ...
                                    agents_ACC(i).iterations(k).info.optimized;
            
            no_cons_used_ACC(k,i) = agents_ACC(i).iterations(k).info.num_cons;
            no_cons_active(k,i) = size(agents_ACC(i).iterations(k).active_deltas, 1);

        end
        p.ping();
    end

end

%% plot
fig = initfig('ACC iterations', 4);
ax = subplot(211, 'YScale', 'log');
grid on
hold on
hs1 = plot(convergence_ACC, '-x', 'color', green);
ylabel('|f(x_k^i) - f(x^*) |')
title(sprintf('ACC convergence for AC OPF %s with %i scenarios', ac.model_name, N));

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility_ACC, '-x', 'color', green);

ylabel('% violated');
xlabel('iterations');

initfig('ACC timing', 2);
ax = subplot(211);
grid on
hold on
hs1 = plot(time_per_iteration_ACC, 'o', 'color', green);
plot(timing_ACC, 'o-', 'color', green);
ylabel('Time per iteration');

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(optimizations_run_ACC, 'color', green);
ylabel('Optimizations run');
xlabel('Iteration');

initfig('No of constraints', 3);
hs1 = plot(no_cons_used_ACC, 'color', green);
hs2 = plot(no_cons_active, 'color', blue);
hs3 = plot(ones(K,1)*length(C_all), '--');
xlabel('Iteration');
legend([hs1(1) hs2(2) hs3], '|L|', '|A|', '|C_{all}|');
figure(fig);
uistack(fig);

%% make assertions
verify(all_close(xstar, xstar_cent, 1e-4), 'Not optimal');
verify(sum(feasibility_ACC(end,:)) == 0, 'Not all feasible');


