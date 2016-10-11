%% Optimize DC OPF with Yalmip
clear all
clc
addpath('../misc');
%% load model and dimensions
p = progress('Loading model', 1);
m = DC_model('case30');
p.update();
N_t = 24;
%% Define optimization variable

P_G = sdpvar(m.N_G, N_t, 'full');         % real power production per gen and time

%% Define objective function

p = progress('Defining objective function', N_t*m.N_G);
Objective = 0;

% loop over time
for t = 1:N_t
    
    % loop over generators
    for k = 1:m.N_G
        Objective = Objective + m.c_qu(k) * (P_G(k, t))^2 + ...
                                m.c_li(k) * P_G(k, t);
        p.update();
    end
end
%% Define constraints
Constraints = [];

% Remove last row and column to get B_bustilde and invert
B_bustilde = m.B_bus(1:end-1, 1:end-1);
B_bustildeinverted = inv(B_bustilde);

p = progress('Defining constraints', N_t);
% loop over time
for t = 1:N_t
    
    % power injection vector 
    P_inj = (m.C_G * P_G(:, t) - m.P_D(:,t));
    
    % power balance constraints
    Constraints = [Constraints, ...
                   ones(1, m.N_b) * P_inj == 0];
    
    % generator limits
    Constraints = [Constraints, ...
                   m.P_Gmin <= P_G(:, t) <= m.P_Gmax];

    % line flow limits
    Constraints = [Constraints, ...
                   - m.P_fmax <=  ...
                   m.B_f * [B_bustildeinverted * P_inj(1:end-1); 0] ...
                   <= m.P_fmax];
               
    p.update();
end

%% Optimize

fprintf('Optimizing ... ')
optimize(Constraints, Objective)
P_Gopt = value(P_G);

%% Plot results
dock
clf
subplot(111)
bar(1:N_t, P_Gopt', 'stack');
xlabel('Time [h]');
ylabel('Real power production [MW]');
title('Real power production over the day-ahead optimization horizon');
grid on