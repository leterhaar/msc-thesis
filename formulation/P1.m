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

figure(1);
ac.draw_network();

N_t = 24;
wind = wind_model(ac, N_t, 0.2);

% define sample complexity
N = 20;
wind.dummy(N);
% wind.use_forecast();
% wind.dummy(N);
%% Define problem
t = 17; % for now, do a loop over 24 hours later

W_f = sdpvar(2*ac.N_b); 
% W_f is a symmetric real valued matrix
W_s = sdpvar(2*ac.N_b, 2*ac.N_b, N);
% W_s is a set of N symmetric real valued matrices

% maximum up and downspinning vectors
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);

% distribution vectors
d_us = sdpvar(ac.N_G, 1);
d_ds = sdpvar(ac.N_G, 1);

% epigraph notation for alpha
alpha = sdpvar(ac.N_G, 1);


figure(2);
wind.plot(t);
pause(0.001); % to show plot
%% Define objective

% P_G for every generator
Obj = sum(alpha);

% Reserve requirements
lambda = 1;
Obj = Obj + lambda*(ac.c_us' * R_us + ac.c_ds' * R_ds);         

%% Define constraints
C = [];
k_ref = ac.N_b + ac.refbus;


% Deterministic constraints
for k = 1:ac.N_b
    % P_inj (1)
    C = [C, ac.P_min(k) ...
        <= trace(ac.Y_k(k)*W_f) + ac.P_D(t, k) - ac.C_w(k)*wind.P_wf(t) <= ...
           ac.P_max(k)];

    % Q_inj (2)
    C = [C, ac.Q_min(k) ...
        <= trace(ac.Ybar_k(k)*W_f) + ac.Q_D(t, k) <=...
           ac.Q_max(k)];

    % V_bus (3)
    C = [C, ac.V_min(k)^2 ...
        <= trace(ac.M_k(k)*W_f) <= ...
           ac.V_max(k)^2];
end

for i = 1:N
        
    for k = 1:ac.N_b
        % P_inj (1)
        C = [C, ac.P_min(k) ...
            <= trace(ac.Y_k(k)*W_s(:,:,i)) + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i) <= ...
               ac.P_max(k)];

        % Q_inj (2)
        C = [C, ac.Q_min(k) ...
             <= trace(ac.Ybar_k(k)*W_s(:,:,i)) + ac.Q_D(t, k) <= ...
                ac.Q_max(k)];

        % V_bus (3)
        C = [C, ac.V_min(k)^2 ...
             <= trace(ac.M_k(k)*W_s(:,:,i)) <= ...
                ac.V_max(k)^2];
    end
    

    for j = 1:ac.N_G

        % bus index
        k = ac.Gens(j);

        % Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= trace(ac.Y_k(k)*(W_s(:,:,i)-W_f)) - ac.C_w(k)*wind.P_m(t, i) <= ...
                R_us(j)];
            
        % relate W_s and W_f through d_ds and d_us
        C = [C, trace(ac.Y_k(k)*(W_s(:,:,i)-W_f)) - ac.C_w(k)*wind.P_m(t, i) == ...
            d_us(j) * max(0, -wind.P_m(t, i)) - d_ds(j) * max(0, wind.P_m(t, i))];
    end
    
    C = [C, W_s(:,:,i) >= 0];
    C = [C, W_s(k_ref, k_ref, i) == 0];

end

% PSD constraint on W_f 
C = [C, W_f >= 0];

% refbus constraint
C = [C, W_f(k_ref, k_ref) == 0];

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];

C = [C, ones(1, ac.N_G)*d_us == 1, ones(1, ac.N_G)*d_ds == 1];

%% Optimize
opt = sdpsettings('verbose', 0, 'debug', 1, 'solver', 'mosek');
diagnostics = optimize(C, Obj, opt);

%% Evaluate
Wf_opt = value(W_f);
Ws_opt = value(W_s);

Rus_opt = value(R_us);
Rds_opt = value(R_ds);

dus_opt = value(d_us);
dds_opt = value(d_ds);

% if only up or downspinning is required, replace the NaNs with zeros
if any(isnan(dus_opt))
    dus_opt = zeros(ac.N_G, 1);
end
if any(isnan(dds_opt))
    dds_opt = zeros(ac.N_G, 1);
end

% replace really small values with zeros
dds_opt(abs(dds_opt) < 1e-4) = 0;
dus_opt(abs(dus_opt) < 1e-4) = 0;

decided_vars = {Wf_opt, Ws_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};

% calculate R for every scenario
R = zeros(ac.N_G, N);
P_G = zeros(ac.N_G, 1);
for i = 1:N
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R(j, i) = trace(ac.Y_k(k) * (Ws_opt(:,:,i)-Wf_opt)) - ac.C_w(k)*wind.P_m(t, i);
        
        % calculate P_G
        if i == 1
            P_G(j) = trace(ac.Y_k(k) * Wf_opt) + ac.P_D(t, k) - ac.C_w(k)*wind.P_wf(k);
        end
    end
end

% output ranks and table with results
format_result(ac, wind, t, decided_vars, R);

if diagnostics.problem ~= 0
    cprintf('red', '%s \t in %g seconds\n\n', diagnostics.info, ... 
                                                   diagnostics.solvertime);
else
    fprintf('%s \t in %g seconds\n\n',diagnostics.info, ...
                                                   diagnostics.solvertime);
end