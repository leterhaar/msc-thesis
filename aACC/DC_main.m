%% DC ACC algorithm

% add helper path
addpath('../misc');
addpath('../formulation');
addpath('../wind');

%% initialize models

% initialize network model
dc = DC_model('case14');
dc.set_WPG_bus(9);

% initialize wind model
wind = wind_model(dc, 24, 0.2);

% 
epsilon = 5e-2;                         % violation parameter
zeta = 5*dc.N_G;                        % Helly-dimsension
beta = 1e-5;                            % confidence parameter

% determine number of scenarios to generate based on Eq (2-4)
N = ceil(2/epsilon*(zeta-1+log(1/beta)));
N = 100;

% generate scenarios
wind.generate(N);
t_wind = 8;

% divide scenarios and initialize agents
N_agents = 5;
assert(N >= N_agents, 'There cannot be more agents than scenarios');
cut_index = ceil(linspace(1, N+1, N_agents+1));

% generate random connection graph with fixed diameter
dm = 2;
G = random_graph(N_agents, dm, 'rand');
% plot(digraph(G))
%% create and init agents
prg = progress('Initializing', N_agents);
for i = 1:N_agents
    agents(i) = DC_agent(dc, wind, t_wind, cut_index(i), cut_index(i+1)-1); 
    prg.ping();
end
%% 
ngc = ones(N_agents, 1);
t = 1;
infeasible = 0;

% while not converged and still feasible
while all(ngc < 2*dm+1) && not(infeasible)
    prg = progress(sprintf('Iteration %i',t), N_agents);
    
    % loop over agents
    for i = 1:N_agents
        
        % loop over all incoming agents to agent i
        for j = find(G(:, i))';
            
            % add incoming A and J
            agents(i).build(agents(j).A{t}, agents(j).J(t), ...
                                                        agents(j).x(:,t));
            
        end
            
        % update agent
        agents(i).update();
        
        if isinf(agents(i).J(t+1))
            infeasible = 1;
            fprintf('Reached infeasibility');
            break
        end
            
        
        % check if J(t+1) is equal to J(t)
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
xstar = zeros(5*dc.N_G, N_agents);
Js = zeros(t, N_agents);
for i = 1:N_agents
    xstar(:, i) = value(agents(i).x_var);
    Js(:, i) = [agents(i).J]';
end

all_agents_are_close = 1;
for i = 1:N_agents-1
    if not( all_close(xstar(:, i), xstar(:, i+1), 1e-4) )
        all_agents_are_close = 0;
    end
end
if all_agents_are_close
    fprintf('\nAll agents are close\n');
else
    fprintf('\n(!) Not all agents are close\n');
end

% check the solution against all constraints a posteriori
tic
x = sdpvar(5*dc.N_G, 1, 'full');

C_all = DC_f_0(x, dc, wind, t_wind);

for i = 1:N
        C_all = [C_all, DC_f_ineq(x, i, dc, wind, t_wind)];
end

feasible_for_all = 1;
for i = 1:N_agents
    assign(x, xstar(:,i));
    if any(check(C_all) < -1e-6)
        feasible_for_all = 0;
    end
end
if feasible_for_all
    fprintf('All solutions are feasible for all original constraints\n')
else
    fprintf('(!) Some solution is infeasible for all original constraints\n');
end

% calculate central solution
Obj = DC_f_obj(x, dc, wind, t_wind);

opt = sdpsettings('verbose', 0);
optimize(C_all, Obj, opt);
xstar_centralized = value(x);
toc

if all_close(xstar(:,1), xstar_centralized, 1e-4)
    fprintf('Decentralized and centralized solution are the same\n');
else
    fprintf('Central solution is different\n');
end
        
%% plot Js
figure(1)
clf
dock

hold on
grid on
xlabel('iteration')
ylabel('J(x*)');
title('Objective vs iterations');
plot(1:t, value(Obj)*ones(1,t), '-.', 'linewidth', 1.2);
plot(Js, 'r-', 'linewidth', 1);
set(gca, 'xtick', 1:t);
legend('Centralized', 'Agents', 'location', 'se');

