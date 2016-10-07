%% Optimize DC OPF with uncertain WPG using Yalmip
clc
addpath('../wind/')
addpath('../misc/')
%% load model and dimensions
p = progress('Loading model', 2);       % "progress" is used for outputting
                                        % progress of loops to command line
                                        % does not do anything else
% number of timesteps
N_t = 1;

% DC network model
m = DC_model('case14');
p.ping();
% attach to bus 22
m.set_WPG_bus(9);

% Wind model
wind = wind_model(m, 24, 0.2);
p.ping();
%% Define decision variables

P_G = sdpvar(m.N_G, N_t, 'full');       % real power production 
R_us = sdpvar(m.N_G, N_t, 'full');      % upspinning reserve requirement
R_ds = sdpvar(m.N_G, N_t, 'full');      % downspinning reserve requirement
d_us = sdpvar(m.N_G, N_t, 'full');      % upspinning distribution vector
d_ds = sdpvar(m.N_G, N_t, 'full');      % downspinning distribution vector

%% Generate scenarios
epsilon = 5e-2;                         % violation parameter
zeta = m.N_G;                           % Helly-dimsension
beta = 1e-5;                            % confidence parameter

% determine number of scenarios to generate based on Eq (2-4)
N = ceil(2/epsilon*(zeta-1+log(1/beta)));
N = 4;
% simulate and generate
p = progress(sprintf('Generating %d wind scenarios', N), 1);
wind.dummy(N);
P_wf = wind.P_wf;
P_wscens = wind.P_w;
p.ping();
%% Define objective function

p = progress('Defining objective function', m.N_G * N_t);
Objective = 0;

% loop over time
for t = 1:N_t
    
    % loop over generators and add generator cost at time t
    for k = 1:m.N_G
        Objective = Objective + m.c_qu(k) * (P_G(k, t))^2 + ...
                                m.c_li(k) * P_G(k, t);
        p.ping();
    end
    
    % add reserve requirements costs
    Objective = Objective + m.c_us' * R_us(:, t) + m.c_ds' * R_ds(:, t);
    
end
%% Define constraints
Constraints = [];

p = progress('Defining deterministic constraints', N_t);
% loop over time
for t = 1:N_t
    
    % define deterministic power injection vector  
    P_injf = (m.C_G * P_G(:, t) - m.P_D + m.C_w * P_wf(t));

    % power balance constraints
    Constraints = [Constraints, ...
                   ones(1, m.N_b) * P_injf == 0];
    
    % generator limits
    Constraints = [Constraints, ...
                   m.P_Gmin <= P_G(:, t) <= m.P_Gmax];

    % line flow limits
    Constraints = [Constraints, ...
                   - m.P_fmax <=  ...
                   m.B_f * [m.B_bustildeinv * P_injf(1:end-1); 0] ...
                   <= m.P_fmax];
               
    % Non-negativity constraints for reserve requirements       
    Constraints = [Constraints, R_us(:, t) >= 0, R_ds(:, t) >= 0];
    
    p.ping();
    
end
%% Define scenario constraints
p = progress('Defining scenario constraints', N*N_t);
% loop over scenarios
for s = 1:N
    P_w = P_wscens(:, s);
    P_m = P_w - P_wf;
  
    % loop over time
    for t = 1:N_t
        
        % define reserve power
        R = d_us(:, t) * max(0, -P_m(t)) - d_ds(:, t) * max(0, P_m(t));
        
        % define scenario power injection vector
        P_injs = m.C_G * (P_G(:, t) + R) + m.C_w * P_w(t) - m.P_D;
        
        % power balance constraints
        Constraints = [Constraints, ...
                   ones(1, m.N_b) * P_injs == 0];
               
        % generator limits
        Constraints = [Constraints, ...
                       m.P_Gmin <= P_G(:, t) + R <= m.P_Gmax];

        % line flow limits
        Constraints = [Constraints, ...
                       - m.P_fmax <=  ...
                       m.B_f * [m.B_bustildeinv * P_injs(1:end-1); 0]...
                       <= m.P_fmax];
                   
        % reserve requirements constraints
        Constraints = [Constraints, -R_ds(:, t) <= R <= R_us(:, t)];
        p.ping();
    end
end
%% Optimize
opt = sdpsettings('gurobi.TimeLimit', 60);
optimize(Constraints, [], opt);
P_Gopt = value(P_G);
%% Plot results
dock
clf
subplot(111)
bar(1:N_t, P_Gopt', 'stack');
xlabel('Time [h]');
ylabel('Real power production [MW]');
title('Real power production over the day-ahead optimization horizon');
labels = cell(m.N_G,1);
for i = 1:m.N_G
    labels{i} = ['Generator ' int2str(i)];
end
grid on
legend(labels, 'location', 'northwest');

%% Test and validate

% retrieve optimized values and initiate simulator
P_G = value(P_G);
d_ds = value(d_ds);
d_us = value(d_us);
DCsim = DC_simulator(m, P_G, d_ds, d_us, N_t, 1e-5);

% generate large amount of wind scenarios
Nsim = 1e4;
[P_wfsim, ~, P_wscenssim] = wind.generate(Nsim);

p = progress('Simulating dynamic reserve scheduling', Nsim);
for i = 1:Nsim
   DCsim.simulate(P_wscenssim(:, i), P_wfsim);
   p.ping();
end
%% Plot probability on constraint violation per hour
figure(2)
dock
clf

subplot(211)
grid on
hold on

title('Generator limits');
plot([28 29], [epsilon*100 epsilon*100], 'r--');
bar(DCsim.V_gen / DCsim.N * 100);
ylabel('Probability of violation (%)');
legend('Theoretical level');
xlim([0.5,24.5]);
set(gca,'xtick',1:24);

subplot(212)
grid on
hold on

title('Line flow limits');
plot([0 25], [epsilon*100 epsilon*100], 'r--');
bar(DCsim.V_pf / DCsim.N * 100);
ylabel('Probability of violation (%)');
xlabel('Time [h]');
xlim([0.5,24.5]);
set(gca,'xtick',1:24);

%% Simulate static reserve strategy

% choose static reserve distribution vector
d_us = zeros(m.N_G, N_t);
d_us(2, :) = 1;
d_ds = d_us;
DCsim_static = DC_simulator(m, P_G, d_ds, d_us, N_t, 1e-5);

p = progress('Simulating static reserve scheduling', Nsim);
for i = 1:Nsim
   DCsim_static.simulate(P_wscenssim(:, i), P_wfsim);
   p.ping();
end
%%
figure(3)
dock
clf

subplot(211)
grid on
hold on

title('Generator limits');
bar(DCsim_static.V_gen / DCsim.N * 100);
ylabel('Probability of violation (%)');
xlim([0.5,24.5]);
ylim([0 100]);
set(gca,'xtick',1:24);

subplot(212)
grid on
hold on

title('Line flow limits');
bar(DCsim_static.V_pf / DCsim.N * 100);
ylabel('Probability of violation (%)');
xlabel('Time [h]');
xlim([0.5,24.5]);
set(gca,'xtick',1:24);