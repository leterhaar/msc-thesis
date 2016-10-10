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
N = 5;
wind.dummy(N);
% wind.use_forecast();
wind.generate(N);
wind2 = copy(wind);
%% Define problem

t = 17; % for now, do a loop over 24 hours later
wind.use_extremes(t);
N = 2;
W_f = sdpvar(2*ac.N_b); 
% W_f is a symmetric real valued matrix
W_m = sdpvar(2*ac.N_b);
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
    Obj = Obj + ac.c_us(j)*trace(ac.Y_k(k)*W_f);
end

% Reserve requirements
lambda = 0.5;
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
    
    W_s = W_f + W_m * wind.P_m(t, i);
        
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
             <= trace(ac.Y_k(k)*(W_m*wind.P_m(t, i))) - ac.C_w(k)*wind.P_m(t, i) <= ...
                R_us(j)];
    end
    
    C = [C, W_s >= 0];
    C = [C, W_s(k_ref, k_ref) == 0];

end

% PSD constraint on W_f 
C = [C, W_f >= 0];

% refbus constraint
C = [C, W_f(k_ref, k_ref) == 0];

Ysum = zeros_like(ac.Y_k(1));
for k = ac.Gens'
    Ysum = Ysum + ac.Y_k(k);
end

% sum Wm = 1
C = [C, trace(Ysum * W_m) == -1];   

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];
%% Optimize
opt = sdpsettings('verbose', 0);
diagnostics = optimize(C, Obj, opt);

%% Evaluate
Wf_opt = value(W_f);
Wm_opt = value(W_m);

Rus_opt = value(R_us);
Rds_opt = value(R_ds); 

% if no W_m is required, replace NaN with zeros
if all(all(isnan(Wm_opt)))
   Wm_opt = zeros(2*ac.N_b);
end

% extract d from W_m by simulating P_m = -1
d = zeros(ac.N_G, 1);
for j = 1:ac.N_G
    k = ac.Gens(j);
    d(j) = trace(ac.Y_k(k)*-Wm_opt) + ac.C_w(k);
    Rds_opt(j) = -trace(ac.Y_k(k)*Wm_opt*max([wind.P_m(t,:) 0]));
    Rus_opt(j) = trace(ac.Y_k(k)*Wm_opt*min([wind.P_m(t,:) 0]));
end

% normalize such that sum(d) = 1
dus_opt = normV(d);
dds_opt = normV(d); % downspinning = upspinning in this formulation

% create scenario W_s for every mismatch
Ws_opt = zeros(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
   Ws_opt(:,:,i) = Wf_opt + Wm_opt * wind.P_m(t,i);
end

decided_vars = {Wf_opt, Ws_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};

% calculate R for every scenario
R = zeros(ac.N_G, N);
% R2 = zeros(ac.N_G, N);
for i = 1:N
    for j = 1:ac.N_G
        k = ac.Gens(j);
%         R2(j, i) = trace(ac.Y_k(k) * (Ws_opt(:,:,i)-Wf_opt)) - ac.C_w(k)*wind.P_m(t, i);
        R(j, i) = trace(ac.Y_k(k) * (Wm_opt * wind.P_m(t,i))) - ac.C_w(k)*wind.P_m(t, i);
    end
end

% output ranks and table with results
format_result(ac, wind, t, decided_vars, R);

if diagnostics.problem ~= 0
    fprintf('%s (!) \t\t in %g seconds\n\n', diagnostics.info, diagnostics.solvertime);
else
    fprintf('%s \t\tin %g seconds\n\n',diagnostics.info, diagnostics.solvertime);
end

decided_vars{2} = Wm_opt;

% formulate a new problem to check against all PSD constraints
p = formulate('P3');
p.evaluate(ac, wind2, t, decided_vars);