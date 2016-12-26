%% DC formulation for IAPG problem
clear
yalmip('clear');
addpath('../networks/');
addpath('../wind');
addpath('../formulation_dc');
addpath('../experiment');
addpath('../misc');

N_t = 24;      % number of time steps
N = 50;       % number of scenarios
max_its = 500; % maximum number of iterations

% initialize models
dc = DC_model('case14a');
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
fs = cell(N,1);
grad_fs = cell(N,1);
for i = 1:N
    C = DC_cons_scen(x, dc, wind.slice(i));
    C_cent = [C_cent, C];
    constraints{i} = C;
    fs{i} = @(x) DC_f(x, dc, wind);
    grad_fs{i} = @(x) DC_gradient_f(x, dc, wind);
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

[~, its_IAPG] = IAPG(x, fs, grad_fs, ...
                             constraints, ...
                             'x0', x0, ...
                             'default_constraint', C_det, ...
                             'opt_settings', opt_settings, ...
                             'max_its', max_its, ...
                             'verbose', 1);
                         
[~, its_IAPG_light] = IAPG_light(x, fs, grad_fs, ...
                             constraints, ...
                             'x0', x0, ...
                             'default_constraint', C_det, ...
                             'opt_settings', opt_settings, ...
                             'max_its', max_its, ...
                             'verbose', 1, ...
                             'b', N);
                                                  

[~, its_IAPG_register] = IAPG_light_register(x, fs, grad_fs, ...
                             constraints, ...
                             'x0', x0, ...
                             'default_constraint', C_det, ...
                             'opt_settings', opt_settings, ...
                             'max_its', max_its, ...
                             'verbose', 1, ...
                             'b', N);                         
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
    assign(x, its_IAPG_light(i).x);
    residuals = check(C_cent);
    its_IAPG_light(i).feas = sum(residuals < -tol)/Ncons * 100;
    
     % check IPG
    assign(x, its_IAPG_register(i).x);
    residuals = check(C_cent);
    its_IAPG_register(i).feas = sum(residuals < -tol)/Ncons * 100;
    
    p.ping();
end
%% Plot
initfig('Objectives', 1);
Obj_opt = DC_f(x_star, dc, wind);
differences_objective_IAPG_register = abs([its_IAPG_register(2:end).f] - Obj_opt);
differences_objective_IAPG = abs([its_IAPG(2:end).f] - Obj_opt);
differences_objective_IAPG_light = abs([its_IAPG_light(2:end).f] - Obj_opt);

ax = subplot(211);
hold off
semilogy(differences_objective_IAPG_register, 'linewidth', 2);
hold on
plot(differences_objective_IAPG_light, 'linewidth', 2);
plot(differences_objective_IAPG, 'linewidth', 2);
grid on
ylabel('|f(x)-f(x*)|');
legend('IAPG-l-r', 'IAPG-l', 'IAPG');
title('DC Objective');
% Plot feasibility percentage

ax2 = subplot(212);
plot([its_IAPG_register.feas]);
hold on
title('Feasibility');
plot([its_IAPG_light.feas]);
plot([its_IAPG.feas]);
grid on;
ylabel('% violated');
linkaxes([ax, ax2], 'x');
% 
% norm_subgrads = arrayfun(@(i) norm([its_IAPG(i).subgrad]), 1:max_its);
% 
% ax2 = subplot(212);
% linkaxes([ax, ax2], 'x');
% semilogy(norm_subgrads);
% h = ylabel('$\| \tilde \nabla h_i(x_{k+1})\|$');
% set(h, 'interpreter','latex');
% xlabel('iteration');


%% plot times
initfig('Times', 3);
plot([its_IAPG_register.time]);
plot([its_IAPG_light.time]);
plot([its_IAPG.time]);
legend('IAPG-l-r', 'IAPG-l', 'IAPG');
ylabel('Time per iteration');
xlabel('Iteration');

%%
klaarrr