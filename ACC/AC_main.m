%% DC ACC algorithm

% add helper path
addpath('../experiment');

%% initialize models

N = 5;
Nm = 2;
m = ceil(N/Nm);

init_experiment(...
    'model_name',           'case14a',   ...
    'model_windbus',        9, ...
    'model_formulation',    'P3',       ...
    'wind_N',               N,         ...
    'wind_Nm',              Nm,         ...
    'wind_dummy',           1);

t_wind = 17;

% generate random connection graph with fixed diameter
dm = 2;
% G = ones(2) - diag(ones(2,1));
G = random_graph(m, dm, 'rand');
% plot(digraph(G))
%% create and init agents
clear agents

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
while all(ngc < 2*dm+1) && not(infeasible) && t < 10
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
xstar_cell = cell(4, m);
xstar = zeros(3*(2*ac.N_b)^2+2*ac.N_G, m);
Js = zeros(t, m);
for i = 1:m
    xstar_cell(:, i) = values_cell(agents(i).x);
    xstar(:,i) = [vec(xstar_cell{1, i}); vec(xstar_cell{2, i}); ...
                  vec(xstar_cell{3, i}); xstar_cell{4, i}];
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


% % extract and compare V vectors 
% 
% % preallocate voltage vector
% Xs = zeros(2*ac.N_b, m);
% 
% % preallocate mismatch voltage vector
% Xms = zeros_like(Xs);
% 
% % loop over agents and extract
% for i = 1:m
%     W = xstar_cell{1,i};
%     Wm = xstar_cell{2,i};
%     Xs(:, i) = sqrt(diag(W)) .* sign(W(:,1));
%     Xms(:, i) = sqrt(diag(Wm)) .* sign(Wm(:,1));
% end
% 
% for i = 1:m-1
%     for j = i+1:m
%        if not( all_close(Xs(:, i), Xs(:, j), 1e-3) )
%         all_agents_are_close = 0;
%         fprintf('|| X_%i - X_%i || = %g\n', i, j,...
%                                         norm(Xs(:,i)-Xs(:,j)));
%        end 
%        if not( all_close(Xms(:, i), Xms(:, j), 1e-3) )
%         all_agents_are_close = 0;
%         fprintf('|| Xm_%i - Xm_%i || = %g\n', i, j,...
%                                         norm(Xs(:,i)-Xs(:,j)));
%        end 
%     end
% end

%%
% check the solution against all constraints a posteriori
tic
x = {   sdpvar(2*ac.N_b), ...       Wf
        sdpvar(2*ac.N_b), ...       Wmus
        sdpvar(2*ac.N_b), ...       Wmds
        sdpvar(2*ac.N_G, 1)}; ...   Rus and Rds

%% check constraints
feasible_for_all = 1;
for agent_id = 1:m
    for i = 1:N
       problem = AC_check({xstar_cell{:,agent_id}}, ac, wind.slice(i), t_wind);
       if problem
           feasible_for_all = 0;
           fprintf('Problem type %i with agent %i for scenario %i\n', problem, agent_id, i);
       end
    end
end

%%
% construct all scenario constraints
C_all = [];
for i = 1:N
        C_all = [C_all, AC_cons_scen(x, ac, wind.slice(i), t_wind)];
end

C_all = [C_all, AC_cons_det(x, ac, wind.slice(1), t_wind)];

% calculate centralized solution
Obj = AC_f(x, ac, wind, t_wind);

opt = sdpsettings('verbose', 0);
optimize(C_all, Obj, opt);

xstar_cell_c = values_cell(x);
xstar_c = [vec(xstar_cell_c{1}); vec(xstar_cell_c{2}); vec(xstar_cell{3}); xstar_cell_c{4}];
W_c = xstar_cell_c{1};
Wmup_c = xstar_cell_c{2};
Wmds_c = xstar_cell_c{2};

%%
% compare all solutions to the centralized solution
agreement = true;
for i = 1:m
    if not(all_close(xstar(:,i), xstar_c, 1e-3))
        fprintf('|| x_%i - x_c || = %g \n', i, norm(xstar_c-xstar(:,i)));
        agreement = false;
    end
end

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

if agreement
    fprintf('Decentralized and centralized solution are the same\n');
else
    fprintf('(!) Centralized solution is different\n');
end

%% plot disagreement and feasibility
initfig('Objectives', 1);
Obj_opt = AC_f(xstar_cell_c, ac, wind, t_wind);
differences = nan(m, t);
for ag = 1:m
    for k = 1:t
        local_obj = agents(ag).J(k);
        differences(ag, k) = abs(Obj_opt - local_obj);
    end
end

% check feasibility
feas = nan(m, t);
for ag = 1:m
    for k = 1:t
        total_violated = 0;
        for i = 1:N
            if ~isempty(agents(ag).x_hist{k})
                [~, res] = AC_f_check([agents(ag).x_hist{k}], i, ac, wind, t_wind);
                total_violated = total_violated + sum(res < -1e-6);
            end
        end
        feas(ag, k) = total_violated / (N_j*N) * 100;
    end
end

ax = subplot(211);
hold off
semilogy(differences', 'linewidth', 2);
hold on
grid on
ylabel('|f(x_i)-f(x*)|');
legend({'Ag1','Ag2','Ag3'});
title('Objective');
% Plot feasibility percentage

ax2 = subplot(212);
hold on
title('Feasibility');
grid on
plot(feas');
legend({'Ag1','Ag2','Ag3','Ag4','Ag5'});
ylabel('% violated');
linkaxes([ax, ax2], 'x');
xlabel('Iteration');

%% plot Js
initfig('Js',1);
hold on
grid on
xlabel('iteration')
ylabel('J(x*)');
title('Objective vs iterations for ACC');
plot(1:t, value(Obj)*ones(1,t), '-.', 'linewidth', 1.2);
plot(Js, 'r-', 'linewidth', 1);
set(gca, 'xtick', 1:t);
legend('Centralized', 'Agents', 'location', 'se');

%% show image of agents constraints
% enable figure
figure(3);
set(gcf, 'Name', 'Constraint exchange');
N_j = 61;
% make all params
scens = repmat(1:N, N_j, 1);
all_params = [reshape(scens, N_j*N, 1) repmat([1:N_j]', N, 1)];
height = N*N_j;

% loop over agents
for agent_id = 1:m
    
    % preallocate image
    the_image = zeros(height, t);
    
    % enable subfigure
    subplot(1,m,agent_id);
    
    % loop over iterations
    for iteration = 1:t
        
        % retrieve set of constraints
        A = agents(agent_id).A{iteration};
        if isempty(A)
            break;
        end
                
        % loop over pixel rows
        for row = 1:height
            
        % see where difference is 0 (identical)
            same_scen = A(:, 1) - all_params(row, 1) == 0;
            same_j = A(:, 2) - all_params(row, 2) == 0;

            % if this is on the same place, we have a match
            if any(same_scen & same_j)
                
                % set pixel to 1
                the_image(row, iteration) = 1;
            end

        end
        
    end

    % plot picture
    imagesc(the_image);
    
    % set labels etc
    xlabel('Iterations');
    ylabel('Scenario');
    ax = gca;
   
    ax.YTick = ceil(N_j/2):N_j:N*N_j;
    ax.YTickLabels = 1:N;
    singletick
    title(sprintf('Agent %i', agent_id));
    
end
