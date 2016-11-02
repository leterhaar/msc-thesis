%% DC ACC algorithm

% add helper path
addpath('../experiment');

%% initialize models

N = 10;
Nm = 3;

init_experiment(...
    'model_name',           'case14',   ...
    'model_formulation',    'P3',       ...
    'wind_N',               N,         ...
    'wind_Nm',              Nm);

t_wind = 8;

% generate random connection graph with fixed diameter
diam = 2;
G = random_graph(m, diam, 'rand');
% plot(digraph(G))
%% create and init agents
prg = progress('Initializing', m);
for i = 1:m
    agents(i) = AC_agent(ac, wind, t_wind, cut(i,1), cut(i,2)); 
    prg.ping();
end
%% 
ngc = ones(m, 1);
t = 1;
infeasible = 0;

% while not converged and still feasible
while all(ngc < 2*diam+1) && not(infeasible)
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
xstar_cell = cell(3, m);
xstar = zeros(2*(2*ac.N_b)^2+2*ac.N_G, m);
Js = zeros(t, m);
for i = 1:m
    xstar_cell(:, i) = values_cell(agents(i).x);
    xstar(:,i) = [vec(xstar_cell{1, i}); vec(xstar_cell{2, i});  ...
                                                xstar_cell{3, i}];
    Js(:, i) = [agents(i).J]';
end

% compare agent solutions
all_agents_are_close = 1;
for j = 1:m-1
    for i = j+1:m
       if not( all_close(xstar(:, i), xstar(:, j), 1e-3) )
        all_agents_are_close = 0;
        fprintf('Biggest diff between x agent %i and %i: \t %g\n', i, j,...
                                        max(abs(xstar(:,i)-xstar(:,j))));
       end 
    end
end

if all_agents_are_close
    fprintf('\nAll agents are close\n');
else
    fprintf('\n(!) Not all agents are close\n');
end

% extract and compare V vectors 

% preallocate voltage vector
Xs = zeros(2*ac.N_b, m);

% preallocate mismatch voltage vector
Xms = zeros_like(Xs);

% loop over agents and extract
for i = 1:m
    W = xstar_cell{1,i};
    Wm = xstar_cell{2,i};
    Xs(:, i) = sqrt(diag(W)) .* sign(W(:,1));
    Xms(:, i) = sqrt(diag(Wm)) .* sign(Wm(:,1));
end

for i = 1:m-1
    for j = i+1:m
       if not( all_close(Xs(:, i), Xs(:, j), 1e-3) )
        all_agents_are_close = 0;
        fprintf('Biggest diff between V agent %i and %i: \t %g\n', i, j,...
                                        max(abs(Xs(:,i)-Xs(:,j))));
       end 
    end
end

%%
% check the solution against all constraints a posteriori
tic
x = {   sdpvar(2*ac.N_b), ...       Wf
        sdpvar(2*ac.N_b), ...       Wm
        sdpvar(2*ac.N_G, 1)}; ...   Rus and Rdssdpvar(5*dc.N_G, 1, 'full');

% construct all scenario constraints
for i = 1:N
        C_all = [C_all, AC_f_ineq(x, i, ac, wind, t_wind)];
end
%%
feasible_for_all = 1;
for i = 1:m
    N_j = 6*ac.N_b + 2*ac.N_G + 1;
    residuals = zeros(N*N_j,1);
    for j = 1:N
        offset = (j-1)*N_j;
        [~, residuals(offset+1:offset+N_j)] = ...
                        AC_f_check(xstar_cell(:, i), i, ac, wind, t_wind);
    end
    if any(residuals < -1e-3)
        feasible_for_all = 0;
        fprintf('Min residual agent %i: \t %g\n', i, min(residuals));
        C_all(residuals < -1e-3)
        fprintf('\n\n');
    end
end
if feasible_for_all
    fprintf('All solutions are feasible for all original constraints\n')
else
    fprintf('(!) Some solution is infeasible for all original constraints\n');
end

% add f_0
C_all = [C_all, AC_f_0(x, ac, wind, t_wind)];

%%
tic
% calculate centralized solution
Obj = AC_f_obj(x, ac, wind, t_wind);

opt = sdpsettings('verbose', 0);
optimize(C_all, Obj, opt);
xstar_cell_c = values_cell(x);
xstar_c = [vec(xstar_cell_c{1}); vec(xstar_cell_c{2}); xstar_cell_c{3}];
toc
%%
% compare all solutions to the centralized solution
agreement = true;
for i = 1:m
    if not(all_close(xstar(:,i), xstar_c, 1e-3))
        fprintf('Central solution is different, max diff: \t %g \n', ...
                        max(abs(xstar_c-xstar(:,i))));
        agreement = false;
    end
end

if agreement
    fprintf('Decentralized and centralized solution are the same\n');
else
    fprintf('(!) Centralized solution is different\n');
end
%%

% plot Js
figure(2)
clf
dock

hold on
grid on
xlabel('iteration')
ylabel('J(x*)');
title('Objective vs iterations for ACC');
plot(1:t, value(Obj)*ones(1,t), '-.', 'linewidth', 1.2);
plot(Js, 'r-', 'linewidth', 1);
set(gca, 'xtick', 1:t);
legend('Centralized', 'Agents', 'location', 'se');

figure(2)
klaarrr