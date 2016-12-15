%% DC ACC algorithm
clear
yalmip('clear');
% add helper path
addpath('../experiment');
N = 100;
m = 5;
% init experiment
init_experiment(...
    'model_name', 'case14a', ... adapted 14 bus network
    'model_formulation', 'DC', ...
    'model_windbus', 9,...
    'wind_N', N, ... scenarios
    'wind_Nm', ceil(N/m), ... scenarios per agent
    'algorithm', 'ACC', ...
    'wind_dummy', 1);
t_wind = 8;
% generate random connection graph with fixed diameter
dm = 2;
G = random_graph(m, dm, 'rand');
% G = ones(2) - diag(ones(2,1));
%% create and init agents
prg = progress('Initializing', m);
clear agents
for i = 1:m
    agents(i) = DC_agent(dc, wind, t_wind, cut(i, 1), cut(i, 2)); 
    prg.ping();
end
%% 
ngc = ones(m, 1);
t = 1;
infeasible = 0;

% while not converged and still feasible
while any(ngc < 2*dm+1) && not(infeasible)
    prg = progress(sprintf('Iteration %i',t), m);
    
    % loop over agents
    for i = 1:m
        
        % loop over all incoming agents to agent i
        for j = find(G(:, i))';
            
            % add incoming A and J
            agents(i).build(agents(j).A{t}, agents(j).J(t));
            
        end
            
        % update agent
        agents(i).update();
        
        if isinf(agents(i).J(t+1))
            infeasible = 1;
            fprintf('Reached infeasibility');
            break
        end
            
        
        % check if J(t+1) is J(t)
        if all_close(agents(i).J(t+1), agents(i).J(t), 1e-9)
            ngc(i) = ngc(i) + 1;
        else
            ngc(i) = 1;
        end
        
        prg.ping();
    end
    
    % update iteration number     
    t = t + 1;
end

%% Validation 

% store value for J and x for all agents
xstar = zeros(5*dc.N_G, m);
Js = zeros(t, m);
for i = 1:m
    xstar(:, i) = value(agents(i).x);
    Js(:, i) = [agents(i).J]';
end

% compare agent solutions
all_agents_are_close = 1;
for j = 1:m-1
    for i = j+1:m
       if not( all_close(xstar(:, i), xstar(:, j), 1e-3) )
        all_agents_are_close = 0;
        fprintf('|| x_%i - x_%i || = %g\n', i, j,...
                                        norm(xstar(:,i)-xstar(:,j)));
       end 
    end
end


% check the solution against all constraints a posteriori

x = sdpvar(5*dc.N_G, 1, 'full');
C_all = [];

for i = 1:N
        C_all = [C_all, DC_f_ineq(x, i, dc, wind, t_wind)];
end

feasible_for_all = 1;
for i = 1:m
    N_j = 4*dc.N_G + 2*dc.N_l;
    residuals = zeros(N*N_j,1);
    for j = 1:N
        offset = (j-1)*N_j;
        [~, residuals(offset+1:offset+N_j)] = ...
                        DC_f_check(xstar(:, i), i, dc, wind, t_wind);
    end
    if any(residuals < -1e-6)
        feasible_for_all = 0;
        fprintf('Min residual agent %i: \t %g\n', i, min(residuals));
        fprintf('\n\n');
    end

        
end

    
%%
% calculate central solution
total_scenario_constraints = length(C_all);
C_all = [C_all, DC_f_0(x, dc, wind, t_wind)];
Obj = DC_f_obj(x, dc, wind, t_wind);

opt = sdpsettings('verbose', 0);
optimize(C_all, Obj, opt);
xstar_centralized = value(x);

central_local_same = 1;
for i = 1:m
   if not( all_close(xstar(:, i), xstar_centralized, 1e-3) )
    central_local_same = 0;
    fprintf('|| x_%i - x_c || =  %g\n', i,...
                                    norm(xstar(:,i)-xstar_centralized));
   end 
end

% print messages
if all_agents_are_close
    fprintf('\nAll agents are close\n');
else
    fprintf('\n(!) Not all agents are close\n');
end
if feasible_for_all
    fprintf('All solutions are feasible for all original constraints\n')
else
    fprintf('(!) Some solution is infeasible for all original constraints\n');
end
if central_local_same
    fprintf('Decentralized and centralized solution are the same\n');
else
    fprintf('(!) Central solution is different\n');
end
        
%% plot disagreement and feasibility
initfig('Objectives', 1);
Obj_opt = DC_f_obj(xstar_centralized, dc, wind, t_wind);
differences = nan(m, t);
for ag = 1:m
    for k = 1:t
        local_obj = DC_f_obj([agents(ag).x_hist(:, k)], dc, wind, t_wind);
        differences(ag, k) = abs(Obj_opt - local_obj);
    end
end

% check feasibility
DC_active_constraints = nan(m,t);
feas = nan(m, t);
for ag = 1:m
    for k = 1:t
        total_violated = 0;
        for i = 1:N
            [~, res] = DC_f_check([agents(ag).x_hist(:,k)], i, dc, wind, t_wind);
            total_violated = total_violated + sum(res < -1e-6);
        end
        feas(ag, k) = total_violated / (N_j*N) * 100;
        DC_active_constraints(ag, k) = size(agents(ag).A{k}, 1);
    end
end

ax = subplot(211);
hold off
semilogy(differences', 'linewidth', 2);
hold on
grid on
ylabel('|f(x_i)-f(x*)|');
% legend({'Ag1','Ag2','Ag3','Ag4','Ag5', 'Ag6', 'Ag7', 'Ag8'});
title('Objective');
% Plot feasibility percentage

ax2 = subplot(212);
hold on
title('Feasibility');
grid on
plot(feas');
% legend({'Ag1','Ag2','Ag3','Ag4','Ag5', 'Ag6', 'Ag7', 'Ag8'});
ylabel('% violated');
linkaxes([ax, ax2], 'x');
xlabel('Iteration');
% 
% norm_subgrads = arrayfun(@(i) norm([its_IAPG(i).subgrad]), 1:max_its);
% 
% ax2 = subplot(212);
% linkaxes([ax, ax2], 'x');
% semilogy(norm_subgrads);
% h = ylabel('$\| \tilde \nabla h_i(x_{k+1})\|$');
% set(h, 'interpreter','latex');
% xlabel('iteration');
%% Same, but then with the new algorithm

% build problem
x_sdp = sdpvar(5*dc.N_G,1, 'full');
delta_sdp = sdpvar(1,3, 'full');
f = @(x) DC_f_obj(x, dc, wind, t_wind);
default_constraint = DC_f_0(x_sdp, dc, wind, t_wind);
constraints_delta = DC_f_ineq_delta(x_sdp, delta_sdp, dc, t_wind);
opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');
%%
% build deltas
deltas = [];
for i = 1:N
    deltas = [deltas; wind.P_w(t_wind, i), ...
                  max(0, -wind.P_m(t_wind, i)), ...
                  max(0, wind.P_m(t_wind, i))];
end
[xstar_acc, agents] = ACC(x_sdp, delta_sdp, deltas, f, constraints_delta, ...
                          'verbose', 1, ...
                          'default_constraint', default_constraint, ...
                          'n_agents', m,...
                          'diameter', dm,...
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
p = progress('Checking constraints', m);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence(k,i) = abs(agents(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(x, agents(i).iterations(k).x)
        feasibility(k,i) = sum(check(C_all) < -1e-5) / N * 100;
        
        % store times
        time_per_iteration(k,i) = agents(i).iterations(k).time;
        
        % store no of constraints
        ACC_active_deltas(k,i) = ...
                            size(agents(i).iterations(k).active_deltas, 1);
    end
    p.ping();
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
plot(time_per_iteration);

initfig('Activeness', 3);
plot([zeros(1, m); DC_active_constraints'], '--');
plot(ACC_active_deltas);
plot(repmat(total_scenario_constraints, xlim), ':');
