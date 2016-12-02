%% DC formulation for IAPG problem
clear
yalmip('clear');
addpath('../networks/');
addpath('../wind');
addpath('../formulation_dc');
addpath('../experiment');

N_t = 24;       % number of time steps
N = 50;        % number of scenarios
max_its = 250; % maximum number of iterations

% initialize models
dc = DC_model('case14a');
dc.set_WPG_bus(9);
wind = wind_model(dc, N_t, 0.2);
wind.generate(N);

opt_settings = sdpsettings('verbose', 0);
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

info = optimize(C_cent, Obj, opt_settings);
assert(not(info.problem), info.info);
x_star = value(x);

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
                             'verbose', 1);
                         
[x_IAPG2, its_IAPG2] = IAPG_light(x, @(x) DC_f(x, dc, wind), ...
                             @(x) DC_gradient_f(x, dc, wind), ...
                             constraints, ...
                             'x0', x0, ...
                             'default_constraint', C_det, ...
                             'opt_settings', opt_settings, ...
                             'max_its', max_its, ...
                             'verbose', 1, ...
                             'b', N);
                                                  
[x_IPG, its_IPG] = IPG(x, @(x) DC_f(x, dc, wind), ...
                             @(x) DC_gradient_f(x, dc, wind), ...
                             constraints, ...
                             'x0', x0, ...
                             'default_constraint', C_det, ...
                             'opt_settings', opt_settings, ...
                             'max_its', max_its, ...
                             'verbose', 1);
                         
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
    assign(x, its_IAPG2(i).x);
    residuals = check(C_cent);
    its_IAPG2(i).feas = sum(residuals < -tol)/Ncons * 100;
    
     % check IPG
    assign(x, its_IPG(i).x);
    residuals = check(C_cent);
    its_IPG(i).feas = sum(residuals < -tol)/Ncons * 100;
    
    p.ping();
end
%% plot objectives
initfig('DC Objectives', 4);
Obj_opt = DC_f(x_star, dc, wind);
differences_objective_IPG = abs([its_IPG(2:end).f] - Obj_opt);
differences_objective_IAPG = abs([its_IAPG(2:end).f] - Obj_opt);
differences_objective_IAPG2 = abs([its_IAPG2(2:end).f] - Obj_opt);
differences_objective_IAPG(1:N) = nan;
% differences_objective_IAPG2(1:N) = nan;
% differences_objective_IPG(1:b) = nan;
hold off
semilogy(differences_objective_IPG, 'linewidth', 2);
hold on
plot(differences_objective_IAPG, 'linewidth', 2);
plot(differences_objective_IAPG2, 'linewidth', 2);
grid on
title('Objectives');
ylabel('|f(x)-f(x*)|');
xlabel('Iteration');
legend('IPG', 'IAPG', 'IAPG light');

%% Plot feasibility percentage
initfig('Feasibility', 5);
plot([its_IPG.feas]);
plot([its_IAPG2.feas]);
plot([its_IAPG.feas]);
legend('IPG', 'IAPG', 'IAPG light');
xlabel('iteration');
ylabel('% violated');

%% plot times
initfig('Times', 6);
plot([its_IPG.time]);
plot([its_IAPG.time]);
plot([its_IAPG2.time]);
legend('IPG', 'IAPG', 'IAPG light');
ylabel('Time per iteration');
xlabel('Iteration'); 