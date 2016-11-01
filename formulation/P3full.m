%% AC OPF
% OtH 6-9-16
init_experiment( ...
    'model_name',       'case14', ...
    'model_formulation','P3full', ...
    'model_windbus',    9);
wind2 = copy(wind);
wind.use_extremes(t);
N = 2;
%% Define problem

t = 17; % for now, do a loop over 24 hours later

W_f = sdpvar(2*ac.N_b); 
% W_f is a symmetric real valued matrix
W_mus = sdpvar(2*ac.N_b);
W_mds = sdpvar(2*ac.N_b);
% W_m is a symmetric real valued matrix

% maximum up and downspinning vectors
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);

figure(2);
wind.plot(t);
pause(0.001); % to show plot

%% Define objective
Obj = 0;

% P_G for every generator
for j = 1:ac.N_G
    k = ac.Gens(j);
    Obj = Obj + ac.c_us(j)*(trace(ac.Y_k(k)*W_f)+ac.P_D(t,k));
end

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
%%
for i = 1:N
    
    W_s = W_f + W_mus * max(0, -wind.P_m(t, i)) ...
              - W_mds * max(0, wind.P_m(t,i));
          
          
            
    for k = 1:ac.N_b
        % P_inj (1)
        C = [C, ac.P_min(k) ...
            <= trace(ac.Y_k(k)*W_s) + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i) <= ...
               ac.P_max(k)];

        % Q_inj (2)
        C = [C, ac.Q_min(k) ...
             <= trace(ac.Ybar_k(k)*W_s) + ac.Q_D(t, k) <= ...
                ac.Q_max(k)];

        % V_bus (3)
        C = [C, ac.V_min(k)^2 ...
             <= trace(ac.M_k(k)*W_s) <= ...
                ac.V_max(k)^2];
    end
    

    for j = 1:ac.N_G

        % bus index
        k = ac.Gens(j);

        % Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*wind.P_m(t, i) <= ...
                R_us(j)];
    end
    
    C = [C, W_s >= 0];
    

end

% PSD constraint on W_f 
C = [C, W_f >= 0];

% refbus constraint
C = [C, W_f(k_ref, k_ref) == 0];
C = [C, W_mus(k_ref, k_ref) == 0];
C = [C, W_mds(k_ref, k_ref) == 0];

Ysum = zeros_like(ac.Y_k(1));
for k = ac.Gens'
    Ysum = Ysum + ac.Y_k(k);
end

% sum Wm = 1
C = [C, trace(Ysum * W_mus) == 1];   
C = [C, trace(Ysum * W_mds) == 1];

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];
%% Optimize
opt = sdpsettings('verbose', 0);
diagnostics = optimize(C, Obj, opt);

%% Evaluate

Wf_opt = value(W_f);
Wmus_opt = zero_for_nan(value(W_mus));
Wmds_opt = zero_for_nan(value(W_mds));

Rus_opt = value(R_us);
Rds_opt = value(R_ds); 

% extract d from W_m by simulating P_m = -1
dus_opt = zeros(ac.N_G, 1);
dds_opt = zeros_like(dus_opt);
for j = 1:ac.N_G
    k = ac.Gens(j);
    dus_opt(j) = trace(ac.Y_k(k)*Wmus_opt) + ac.C_w(k);
    dds_opt(j) = trace(ac.Y_k(k)*Wmds_opt) + ac.C_w(k);
end

% create scenario W_s and R for every scenario
Ws_opt = zeros(2*ac.N_b, 2*ac.N_b, N);
R = zeros(ac.N_G, N);
for i = 1:N
    Ws_opt(:,:,i) = Wf_opt + Wmus_opt * max(0, -wind.P_m(t, i)) ...
                          - Wmds_opt * max(0, wind.P_m(t,i));
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R(j, i) = trace(ac.Y_k(k) * (Ws_opt(:,:,i)-Wf_opt)) ...
                                        - ac.C_w(k)*wind.P_m(t, i);
    end
end

decided_vars = {Wf_opt, Ws_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};

% output ranks and table with results
format_result(ac, wind, t, decided_vars, R);

if diagnostics.problem ~= 0
    fprintf('%s (!) \t\t in %g seconds\n\n',  ...
                            diagnostics.info, diagnostics.solvertime);
else
    fprintf('%s \t\tin %g seconds\n\n', ...
                            diagnostics.info, diagnostics.solvertime);
end
