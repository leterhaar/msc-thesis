%% DC formulation for IAPG problem
clear
yalmip('clear');
addpath('../networks/');
addpath('../wind');
addpath('../formulation_dc');
addpath('../experiment');
addpath('../misc');

N_t = 24;      % number of time steps
N = 10;       % number of scenarios
max_its = 100; % maximum number of iterations
batch_size = 3;

% initialize models
dc = DC_model('case9');
dc.set_WPG_bus(9);
wind = wind_model(dc, N_t, 0.2);
wind.dummy(N);

opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');
%% define problem
x = sdpvar(5*dc.N_G, N_t, 'full');
Obj = DC_f(x, dc, wind);

C_det = DC_cons_det(x, dc, wind);
C_cent = [C_det];
constraints = cell(N,1);
for i = 1:N
    C = DC_cons_scen(x, dc, wind.slice(i));
    C_cent = [C_cent, C];
    constraints{i} = C;
end

% solve centralized problem
tic
info = optimize(C_cent, Obj, opt_settings);
assert(not(info.problem), info.info);
x_star = value(x);
toc

% solve problem w/o objective
info = optimize(C_cent, [], opt_settings);
assert(not(info.problem), info.info);
x0 = value(x);

%% run algorithms

[x_IAPG, its_IAPG] = IAPG(x, @(x) DC_f(x, dc, wind), ...
                             @(x) DC_gradient_f(x, dc, wind), ...
                             constraints, ...
                             'x0', x0, ...
                             'default_constraint', C_det, ...
                             'opt_settings', opt_settings, ...
                             'max_its', max_its, ...
                             'verbose', 1, ...
                             'batch', batch_size);
                         
[x_BIAPG, its_BIAPG] = BIAPG(x, @(x) DC_f(x, dc, wind), ...
                             @(x) DC_gradient_f(x, dc, wind), ...
                             constraints, ...
                             'x0', x0, ...
                             'default_constraint', C_det, ...
                             'opt_settings', opt_settings, ...
                             'max_its', max_its, ...
                             'verbose', 1, ...
                             'batch', batch_size);
%% check feasibility a posteriori
tol = 1e-6;
p = progress('Checking constraints', max_its);
Ncons = length(C_cent);
for i = 1:max_its
    
    % check IAPG
    assign(x, its_IAPG(i).x);
    residuals = check(C_cent);
    its_IAPG(i).feas = sum(residuals < -tol)/Ncons * 100;
    
    % check IAPG2
    assign(x, its_BIAPG(i).x);
    residuals = check(C_cent);
    its_BIAPG(i).feas = sum(residuals < -tol)/Ncons * 100;
    
    p.ping();
end
%% plot objectives
initfig('Batch DC Objectives', 7);
Obj_opt = DC_f(x_star, dc, wind);

% calculate absolute differences
differences_objective_IAPG = abs([its_IAPG(2:end).f] - Obj_opt);
differences_objective_BIAPG = abs([its_BIAPG(2:end).f] - Obj_opt);

hold off
semilogy(differences_objective_IAPG, 'linewidth', 2);
hold on
plot(differences_objective_BIAPG, 'linewidth', 2);
grid on
title('Objectives');
ylabel('|f(x)-f(x*)|');
xlabel('Iteration');
legend('IAPG with larger set', 'B-IAPG');
%% Plot feasibility percentage
initfig('Batch DC Feasibility', 8);
plot([its_IAPG.feas]);
plot([its_BIAPG.feas]);
legend('IAPG with larger set', 'B-IAPG');
xlabel('iteration');
ylabel('% violated');

%% plot times
initfig('Batch DC Times', 9);
plot([its_IAPG.time]);
plot([its_BIAPG.time]);
legend('IAPG with larger set', 'B-IAPG');
ylabel('Time per iteration');
xlabel('Iteration'); 
%%
klaarrr