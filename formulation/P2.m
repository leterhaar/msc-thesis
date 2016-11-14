%% AC OPF
% OtH 6-9-16
addpath('../misc');
addpath('../wind');

figure(1);
set(1, 'name', 'Network');
dock

figure(2);
set(2, 'name', 'Wind');
dock
%% Load models
ac = AC_model('case14');
ac.set_WPG_bus(9);
ac.c_us(3) = 5;

figure(1);
ac.draw_network();

N_t = 24;
wind = wind_model(ac, N_t, 0.2);

% define sample complexity
N = 2; 
% wind.generate(N);
% wind.use_forecast();
wind.dummy(N);
%% Define problem
t = 8; % for now, do a loop over 24 hours later

W_f = sdpvar(2*ac.N_b); 
% W_f is a symmetric real valued matrix
W_m = sdpvar(2*ac.N_b);

% maximum up and downspinning vectors
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);

% distribution vectors
d_us = sdpvar(ac.N_G, 1);
d_ds = sdpvar(ac.N_G, 1);

figure(2);
wind.plot(t);
pause(0.001); % to show plot
%% Define objective
Obj = 0;

% P_G for every generator
for j = 1:ac.N_G
    k = ac.Gens(j);
    Obj = Obj + ac.c_us(j)*trace(ac.Y_k(k)*W_f);
end

% Reserve requirements
lambda = 1;
Obj = Obj + lambda*(ac.c_us' * R_us + ac.c_us' * R_ds);         

%% Define constraints
C = [];
k_ref = ac.N_b + ac.refbus;


% Deterministic constraints
for k = 1:ac.N_b
    % P_inj (1)
    C = [C, ac.P_min(k) ...
        <= trace(ac.Y_k(k)*W_f) + ac.P_D(t,k) - ac.C_w(k)*wind.P_wf(t) <= ...
           ac.P_max(k)];

    % Q_inj (2)
    C = [C, ac.Q_min(k) ...
        <= trace(ac.Ybar_k(k)*W_f) + ac.Q_D(t,k) <=...
           ac.Q_max(k)];

    % V_bus (3)
    C = [C, ac.V_min(k)^2 ...
        <= trace(ac.M_k(k)*W_f) <= ...
           ac.V_max(k)^2];
end

for i = 1:N
        
    for k = 1:ac.N_b
        % P_inj (1')
        C = [C, ac.P_min(k) ...
            <= trace(ac.Y_k(k)*W_f) + ac.P_D(t,k) - ac.C_w(k)*wind.P_wf(t) ...
               + ac.C_G(k, :) * d_us * max(0, -wind.P_m(t, i))  ...
               - ac.C_G(k, :) * d_ds * max(0, wind.P_m(t, i)) <= ac.P_max(k)];
           
        % Q_inj (2')
        C = [C, ac.Q_min(k) ...
            <=  trace(ac.Ybar_k(k) * (W_f + wind.P_m(t,i) * W_m)) + ac.Q_D(t,k) ...
            <= ac.Q_max(k)];
        
        % V_bus (3')
        C = [C, ac.V_min(k)^2 ... 
            <= trace(ac.M_k(k)*(W_f + wind.P_m(t,i) * W_m)) ...
            <= ac.V_max(k)^2]; 
        
    end
    
    

    for j = 1:ac.N_G

        % Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= d_us(j) * max(0, -wind.P_m(t, i)) ...
              - d_ds(j) * max(0, wind.P_m(t, i)) <= ...
                R_us(j)];
            
    end
    
    C = [C, (W_f + wind.P_m(t,i)*W_m) >= 0];
end

% PSD constraint on W_f 
C = [C, W_f >= 0];

% refbus constraint
C = [C, W_f(k_ref, k_ref) == 0];
C = [C, W_m(k_ref, k_ref) == 0];

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];

% d_ds and d_us should always sum to 1
C = [C, ones(1, ac.N_G)*d_us == 1, ones(1, ac.N_G)*d_ds == 1];
Ysum = zeros_like(ac.Y_k(1));
for k = ac.Gens'
    Ysum = Ysum + ac.Y_k(k);
end

% sum Wm = 1
C = [C, trace(Ysum * W_m) == -1];

%% Optimize
opt = sdpsettings('verbose', 0, 'debug', 1, 'solver', 'mosek');
diagnostics = optimize(C, Obj, opt);

%% Evaluate
Wf_opt = value(W_f);
Wm_opt = value(W_m);

Rus_opt = value(R_us);
Rds_opt = value(R_ds);

dus_opt = value(d_us);
dds_opt = value(d_ds);

% calculate R for every scenario
R1 = zeros(ac.N_G, N);
R2 = zeros_like(R1);
P_G = zeros(ac.N_G, 1);
Ws_opt = zeros(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    Ws_opt(:,:,i) = Wf_opt + Wm_opt * wind.P_m(t,i);
    for j = 1:ac.N_G
        k = ac.Gens(j);
        
        % R1 based on W_m
        R1(j, i) = trace(ac.Y_k(k) * (Wm_opt * wind.P_m(t,i))) ...
                            - ac.C_w(k)*wind.P_m(t, i);
        
        % R2 based on d_us, d_ds
        R2(j, i) = dus_opt(j)*max(0, -wind.P_m(t, i)) ...
                   - dds_opt(j)*max(0, wind.P_m(t, i));
        
        % calculate P_G
        if i == 1
            P_G(j) = trace(ac.Y_k(k) * Wf_opt) + ac.P_D(t,k) ...
                            - ac.C_w(k)*wind.P_wf(t);
        end
    end
end
% output ranks and table with results
decided_vars = {Wf_opt, Ws_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};
fprintf('===============\nR1 : base R on W_m\n===============\n');
format_result(ac, wind, t, decided_vars, R1);
fprintf('===============\nR2 : based R on d_us, d_ds\n===============\n');
format_result(ac, wind, t, decided_vars, R2);

if all_close(R1, R2)
    fprintf('\nR1 and R2 are the same\n\n');
else
    fprintf('\nR1 and R2 are different\n\n');
end

if diagnostics.problem ~= 0
    cprintf('red', '%s\n\n', diagnostics.info);
else
    fprintf('%s\n\n',diagnostics.info);
end