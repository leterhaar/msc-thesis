%% DC Consensus algorithm
% Simple consensus algorithm for distributed solving of DC network

% Add paths
addpath('../formulation/');
addpath('../misc');
addpath('../wind');

%% Load models

% load network model
dc = DC_model('case_ieee30');
dc.set_WPG_bus(9);

% load wind model
N_t = 24;           % number of timesteps
penetration = 0.2;  % wind power penetration
all_winds = wind_model(dc, N_t, penetration);

%% Generate scenarios
epsilon = 5e-2;                         % violation parameter
zeta = 5*dc.N_G;                        % Helly-dimsension
beta = 1e-5;                            % confidence parameter

% determine number of scenarios to generate based on Eq (2-4)
N = ceil(2/epsilon*(zeta-1+log(1/beta)));
all_winds.generate(N);

t = 8;              % loop through control horizon later 

% generate wind and 
% N = 40;
% all_winds.dummy(N);

%% Initialize agents
SPA = 100;                   % number of scenarios per agent
N_a = ceil(N/SPA);

prg = progress(sprintf('Initializing agents \t'), N_a);
agents = cell(N_a,1);

agent_wind = struct('P_wf', all_winds.P_wf,...
                    'P_w', zeros(N_t, SPA),...
                    'P_m', zeros(N_t, SPA));

for i = 1:N_a-1
    
    
    agent_wind.P_w = all_winds.P_w(:, SPA*(i-1)+1:SPA*i);
    agent_wind.P_m = agent_wind.P_w - agent_wind.P_wf*ones(1,SPA);
    
    agents{i} = DC_agent(dc, agent_wind, t);
    prg.ping();
end

% last one, until the end
agent_wind.P_w = all_winds.P_w(:, SPA*i+1:end);
agent_wind.P_m = agent_wind.P_w - agent_wind.P_wf*ones(1,size(agent_wind.P_w,2));
agents{i+1} = DC_agent(dc, agent_wind, t);
prg.ping();


%% Update agents
k = 1;

% initialize global variables
Z = zeros(5*dc.N_G, 1);

% initialize array of local variables
X = zeros(5*dc.N_G, N_a);

% stopping criterion: if local and global solution are close for all
% agents, stop the algorithm
eps = 1e-1;
max_iter = 50;

while 1
    prg = progress(sprintf('Iteration %3i \t\t\t', k), N_a);
    
    % average results to get global variable (i.e. all equal weight 1/N)
    Z = mean(X, 2);
        
    % for the first iteration, don't take consensus into account
    if k == 1
        c = 1e99;
    else
        c = 2/(2*(k+1)^2);
    end
    
    all_agents_agree = 1;
    
    % loop over agents
    for i = 1:N_a
        
        % check if agents previous local and global solution agree
        if k > 1
            disagreement = norm(X(:, i) - Z);
            if disagreement > eps
                all_agents_agree = 0;
            end
        else
            % due to initialization, all agents will 'agree' in the first
            % iteration, this is a workaround
            all_agents_agree = 0;
        end

        % update agents local estimate
        X(:, i) = agents{i}.update(Z, c);
        
        % output progress
        prg.ping();
    end
    
    if all_agents_agree || k >= max_iter
        break
    end
    
    k = k + 1;
end
fprintf('\r');
%% Calculate central solution
prg = progress('Calculating central solution', 1);
central = DC_agent(dc, all_winds, t);
central.update(Z, 1e99);
prg.ping();

%% save
save('results/DC_scenarios_case30.mat', 'agents', 'central');