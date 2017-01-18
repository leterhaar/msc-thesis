% initialize problem
% script to test whether the AC OPF with only extremes yields the same
% solution as the original fully defined problem

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
N_t = 24;   % optimization horizon
N = 15;      % number of scenarios used for optimization
t = 5; % timestep used for this demonstration (todo: add for everything)
tol = 5e-4;

% load network and wind models
ac = AC_model('case14');
ac.set_WPG_bus(9);
wind = wind_model(ac, N_t, 0.2);

% generate a number of scenarios
wind.dummy(N);
wind_all = copy(wind);
wind.use_extremes(t);
N = 2;

% optimization settings
ops = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug', 1);

% connectivity matrix
G = ones(N)-diag(ones(N,1));
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
cons_fcn = @(x, delta, j) AC_cons_scen(x, ac, ... fake wind object using delta
                    struct('P_w', ones(N_t,1)*delta(1), 'P_m', ones(N_t,1)*delta(2)), t, j);

% define residual function
residuals = @(x, delta, j_des) -AC_g(x, ac, ... again, fake wind object for delta
                    struct('P_w', ones(N_t,1)*delta(1), 'P_m', ones(N_t,1)*delta(2)), t, j_des);


%% solve centralized problem
C_all = default_constraints;
C_all_delta = default_constraints;
for i = 1:N
    C_all = [C_all, AC_cons_scen(x_cell, ac, wind.slice(i), t)];
    
    C_all_delta = [C_all_delta, cons_fcn(x_cell, deltas(i, :), 0)];
end
status = optimize(C_all, objective_fcn(x_cell), sdpsettings('verbose', 0));
verify(not(status.problem), status.info);
xstar_cent = values_cell(x_cell);

% status = optimize(C_all_delta, objective_fcn(x_cell));
% verify(not(status.problem), status.info);
% xstar_cent_using_delta = values_cell(x_cell);
% 
% verify(all_close(xstar_cent, xstar_cent_using_delta), ...
%                                             'Cons fcn is not equivalent');
%% run algorithm

[xstar, agents_ACCA] = ACCA_fcn_cell(x_cell, deltas, objective_fcn, cons_fcn, ...
                             'default_constraint', default_constraints,...
                             'residuals', residuals,...
                             'stepsize', @(k) 1e2,...
                             'use_selector', true,...
                             'opt_settings', ops,...
                             'tolerance', tol,...
                             'connectivity', G,...
                             'diameter', diam,...
                             'max_its', 50,...
                             'n_agents', N,...
                             'verbose', 1,...
                             'debug', 1,...
                             'x0', []);
%%                         
verify(AC_check(xstar, ac, wind, t) == 0, 'Infeasible solution');
%% calculate convergence and feasibility

%TODO check a posteriori 
K = length(agents_ACCA(1).iterations);
convergence_ACCA = nan(K,N);
convergenceX = nan(K,N);
feasibility_ACCA = nan(K,N);
time_per_iteration_ACCA = nan(K,N);
optimal_objective = objective_fcn(xstar_cent);
optimizations_run_ACCA = zeros(K,1);
no_cons_used_ACCA = nan(K,N);
p = progress('Checking constraints', N);
timing_ACCA = zeros(K,1);

for i = 1:N
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACCA(k,i) = abs(agents_ACCA(i).iterations(k).J ...
                                                - optimal_objective);
        [~, convergenceX(k,i)] = all_close(agents_ACCA(i).iterations(k).x, xstar_cent);
        
        % calculate feasibility percentage
        assign_cell(x_cell, agents_ACCA(i).iterations(k).x)
        feasibility_ACCA(k,i) = sum(check(C_all) < -tol) / length(C_all) * 100;
        
        % store times
        time_per_iteration_ACCA(k,i) = agents_ACCA(i).iterations(k).time;
        timing_ACCA(k) = timing_ACCA(k) + agents_ACCA(i).iterations(k).time;

        
        % store total number of iterations run
        if k > 1
            optimizations_run_ACCA(k) = optimizations_run_ACCA(k) + ...
                                    agents_ACCA(i).iterations(k).info.optimized;
            
            no_cons_used_ACCA(k,i) = agents_ACCA(i).iterations(k).info.num_cons;

        end
    end
    p.ping();
end

%% plot
fig = initfig('ACCA iterations', 1);
ax = subplot(211, 'YScale', 'log');
grid on
hold on
hs1 = plot(convergence_ACCA);
ylabel('|f(x_k^i) - f(x^*) |')
title(sprintf('ACC convergence for AC OPF %s with %i scenarios', ac.model_name, N));
legend('Agent 1', 'Agent 2', 'location', 'ne');
singletick

% ax3 = subplot(212);
% semilogy(convergenceX);
% hold on;
% grid on;
% singletick

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility_ACCA);
singletick
ylabel('% violated');
xlabel('Iteration');

initfig('ACC timing', 2);
ax = subplot(211);
grid on
hold on
hs1 = plot(time_per_iteration_ACCA, 'o');
hs2 = plot(timing_ACCA, 'o-');
legend([hs2; hs1], 'Total', 'Agent 1', 'Agent 2', 'location', 'ne')
ylabel('Time per iteration');
title('Timing');
singletick

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(optimizations_run_ACCA);
title('Optimizations run');
ylabel('Number of optimizations run');
xlabel('Iteration');
singletick

initfig('Comparison of constraints global vs ACC', 3);
hs1 = plot(no_cons_used_ACCA);
hs2 = plot([2:K], ones(K-1,1)*(length(C_all) - length(default_constraints)));
legend([hs2; hs1], 'Total', 'Agent 1', 'Agent 2', 'location', 'southeast')
ylims = ylim;
ylim([0 ylims(2)]);
xlabel('Iteration');
ylabel('Number of constraint used');
singletick

figure(fig);
%% make assertions
verify(all_close(xstar, xstar_cent, 1e-3), 'Not optimal');
verify(sum(feasibility_ACCA(end,:)) == 0, 'Not all feasible');
