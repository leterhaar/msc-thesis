%% DC Consensus algorithm
% Simple consensus algorithm for distributed solving of DC network

% Add paths
addpath('../formulation/');
addpath('../misc');
addpath('../wind');

%% Load models

% load network model
ac = AC_model('case14');
ac.set_WPG_bus(9);

% load wind model
N_t = 24;           % number of timesteps
penetration = 0.2;  % wind power penetration
all_winds = wind_model(ac, N_t, penetration);

%% Generate scenarios
epsilon = 5e-2;                         % violation parameter
zeta = 5*ac.N_G;                        % Helly-dimsension
beta = 1e-5;                            % confidence parameter

% determine number of scenarios to generate based on Eq (2-4)
N = ceil(2/epsilon*(zeta-1+log(1/beta)));
all_winds.generate(N);

t = 7;              % loop through control horizon later 

% generate wind and 
N = 1000;
all_winds.generate(N);
figure(2)
all_winds.plot(t);

%% Initialize agents
SPA = 100;                   % number of scenarios per agent
N_a = ceil(N/SPA);

prg = progress(sprintf('Initializing agents \t'), N_a);
clear agents
agent_wind = struct('P_wf', all_winds.P_wf,...
                    'P_w', zeros(N_t, SPA),...
                    'P_m', zeros(N_t, SPA));

for i = 1:N_a-1
    
    
    agent_wind.P_w = all_winds.P_w(:, SPA*(i-1)+1:SPA*i);
    agent_wind.P_m = agent_wind.P_w - agent_wind.P_wf*ones(1,SPA);
    
    agents(i) = AC_agent(ac, agent_wind, t);
    prg.ping();
end

% last one, until the end
agent_wind.P_w = all_winds.P_w(:, SPA*i+1:end);
agent_wind.P_m = agent_wind.P_w - agent_wind.P_wf*ones(1,size(agent_wind.P_w,2));
agents(i+1) = AC_agent(ac, agent_wind, t);
prg.ping();


%% Update agents
k = 1;

% initialize collection of local variables
X = zeros(8*ac.N_b^2+2*ac.N_G, N_a);

% initialize global variables
Wf_g = zeros(2*ac.N_b);
Wm_g = zeros(2*ac.N_b);
Rus_g = zeros(ac.N_G, 1);
Rds_g = zeros(ac.N_G, 1);
Z = {Wf_g, Wm_g, Rus_g, Rds_g};

% stopping criterion: if local and global solution are close for all
% agents, stop the algorithm
eps = 1e-2;
max_iter = 50;
disagreements = zeros(N_a,1);

while 1
    prg = progress(sprintf('Iteration %3i \t\t\t', k), N_a);
    
    % average results to get global variable (i.e. all equal weight 1/N)
    Zvec = mean(X, 2);
    
    % reconstruct Z as cell from the vectorized Z
    idx = 0;
    for j = 1:4
        sliced = Zvec(idx+1:idx+numel(Z{j}));
        idx = idx + numel(Z{j});
        if size(Z{j}, 2) > 1
            % reconstruct matrix
            Z{j} = mat(sliced);
        else
            Z{j} = sliced;
        end 
    end
        
    % for the first iteration, don't take consensus into account
    if k == 1
        c = 1e99;
    else
        c = 4/(2*(k+1)^2);
    end
    
    all_agents_agree = 1;
    
    % loop over agents
    for i = 1:N_a
        
        % check if agents previous local and global solution agree
        if k > 1
            disagreements(i) = norm(X(:, i) - Zvec);
            if disagreements(i) > eps
                all_agents_agree = 0;
            end
        else
            % due to initialization, all agents will 'agree' in the first
            % iteration, this is a workaround
            all_agents_agree = 0;
        end

        % update agents local estimate
        X(:, i) = agents(i).update(Z, c);
        
        % output progress
        prg.ping();
    end
    
    fprintf('\b max disagreement: %g\r', max(disagreements));
    
    if all_agents_agree || k >= max_iter
        break
    end
    
    k = k + 1;
end
fprintf('\r');
%% Calculate central solution
prg = progress('Calculating central solution', 1);
central = AC_agent(ac, all_winds, t);
central.update(Z, 1e99);
prg.ping();

%% save
save('results/AC_dummy.mat', 'agents', 'central', 'ac');