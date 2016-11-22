%% DC ACC algorithm

% add helper path
addpath('../experiment');
N = 10;
m = 3;
% init experiment
init_experiment(...
    'model_name', 'case14a', ... adapted 14 bus network
    'model_formulation', 'DC', ...
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
while all(ngc < 2*dm+1) && not(infeasible)
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
%         assign(x, xstar(:, i));
%         check(C_all(residuals < -1e-6));
        fprintf('\n\n');
    end
end
%%
% calculate central solution
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
        
        
%% plot Js
initfig('Js', 1);

hold on
grid on
xlabel('iteration')
ylabel('J(x*)');
title('Objective vs iterations for ACC');
plot(1:t, value(Obj)*ones(1,t), '-.', 'linewidth', 2, 'color', green);
plot(Js, '-', 'linewidth', 2, 'color', blue);
set(gca, 'xtick', 1:t);
legend('Centralized', 'Agents', 'location', 'se');


%% show image of agents constraints
% enable figure
figure(2);
set(gcf, 'Name', 'Constraint exchange');

% make all params
scens = repmat(1:N, N_j, 1);
all_params = [reshape(scens, N_j*N, 1) repmat([1:N_j]', N, 1)];
height = N*N_j;

% loop over agents
for agent_id = 1:m
    
    % preallocate image
    image = zeros(height, t);
    
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
                image(row, iteration) = 1;
            end

        end
        
    end

    % plot picture
    imagesc(image);
    
    % set labels etc
    xlabel('Iterations');
    ylabel('Scenario');
    ax = gca;
   
    ax.YTick = ceil(N_j/2):N_j:N*N_j;
    ax.YTickLabels = 1:N;
    singletick
    title(sprintf('Agent %i', agent_id));
    
end

%% plot disagreement
initfig('disagreements', 3);
disagreements = zeros(t, m);

% loop over agents
for i = 1:m
    
    % loop over iterations
    for it = 1:t
        
        disagreements(it, i) = norm(agents(i).x_hist(:, it) - xstar_centralized);
    end
end

hold off
semilogy(disagreements, 'linewidth', 2, 'color', blue);
grid on