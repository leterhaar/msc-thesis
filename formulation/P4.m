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
ac.P_D = ac.P_D .* 2;
N_t = 24;
wind = wind_model(ac, N_t, 0.2);

% define sample complexity 
N = 10;
% wind.dummy(N);
% wind.use_forecast();
wind.generate(N);
% wind2 = copy(wind);
% wind.use_extremes(t);
% N = 2;
%% Define problem

t = 17; % for now, do a loop over 24 hours later

W_f = sdpvar(2*ac.N_b); 
% W_f is a symmetric real valued matrix
W_mus = sdpvar(2*ac.N_b);
W_mds = sdpvar(2*ac.N_b);
% W_m is a symmetric real valued matrix

figure(2);
wind.plot(t);
pause(0.001); % to show plot

%% Define objective
Obj = 0;

% P_G for every generator
lambda = 0.5;
P_mmax = max([wind.P_m(t,:) 0]);
P_mmin = min([wind.P_m(t,:) 0]);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Obj = Obj + ac.c_us(j)*(trace(ac.Y_k(k)*W_f) + lambda*( ...
                            trace(ac.Y_k(k)*W_mds*P_mmax) + ...
                            trace(ac.Y_k(k)*W_mus*P_mmin)));
end
%% Define constraints
C = [];
k_ref = ac.N_b + ac.refbus;

for k = 1:ac.N_b
    % P_inj (1)
    C = [C, ac.P_min(k) ...
        <= trace(ac.Y_k(k)*W_f) + ac.P_D(t,k) - ac.C_w(k)*wind.P_wf(t) <= ...
           ac.P_max(k)];
    % P_inj upspsinning (1')
    C = [C, ac.P_min(k) ...
        <= trace(ac.Y_k(k)*(W_f+W_mus*P_mmin)) + ac.P_D(t,k) - ac.C_w(k)*(P_mmin+wind.P_wf(t)) <= ...
           ac.P_max(k)];
    % P_inj downspinning (1')
    C = [C, ac.P_min(k) ...
        <= trace(ac.Y_k(k)*(W_f+W_mds*P_mmax)) + ac.P_D(t,k) - ac.C_w(k)*(P_mmax+wind.P_wf(t)) <= ...
           ac.P_max(k)];
    

    % Q_inj (2)
    C = [C, ac.Q_min(k) ...
        <= trace(ac.Ybar_k(k)*W_f) + ac.Q_D(t, k) <=...
           ac.Q_max(k)];
    % Q_inj upspinning (2')
    C = [C, ac.Q_min(k) ...
        <= trace(ac.Ybar_k(k)*(W_f+W_mus*P_mmin)) + ac.Q_D(t, k) <=...
           ac.Q_max(k)];
    % Q_inj downspinning (2')
    C = [C, ac.Q_min(k) ...
        <= trace(ac.Ybar_k(k)*(W_f+W_mds*P_mmax)) + ac.Q_D(t, k) <=...
           ac.Q_max(k)];

    % V_bus (3)
    C = [C, ac.V_min(k)^2 ...
        <= trace(ac.M_k(k)*W_f) <= ...
           ac.V_max(k)^2];
    % V_bus upspinning (3')
    C = [C, ac.V_min(k)^2 ...
        <= trace(ac.M_k(k)*(W_f+W_mus*P_mmin)) <= ...
           ac.V_max(k)^2];
    % V_bus downspinning (3)
    C = [C, ac.V_min(k)^2 ...
        <= trace(ac.M_k(k)*(W_f+W_mds*P_mmax)) <= ...
           ac.V_max(k)^2];
end


% PSD constraints
C = [C, W_f >= 0];
C = [C, (W_f + W_mus*P_mmin) >= 0];
C = [C, (W_f + W_mds*P_mmax) >= 0];

% refbus constraints
C = [C, W_f(k_ref, k_ref) == 0];
C = [C, W_mus(k_ref, k_ref) == 0];
C = [C, W_mds(k_ref, k_ref) == 0];
%% Optimize
opt = sdpsettings('verbose', 1);
diagnostics = optimize(C, Obj, opt);

%% Evaluate
Wf_opt = value(W_f);
Wmus_opt = value(W_mus);
Wmds_opt = value(W_mds);

% extract d from W_m by simulating P_m = -1
d_us = zeros(ac.N_G, 1);
d_ds = zeros(ac.N_G, 1);
Rds_opt = zeros(ac.N_G, 1);
Rus_opt = zeros(ac.N_G, 1);

for j = 1:ac.N_G
    k = ac.Gens(j);
    d_us(j) = trace(ac.Y_k(k)*-Wmus_opt) + ac.C_w(k);
    d_ds(j) = trace(ac.Y_k(k)*-Wmds_opt) + ac.C_w(k);
    Rus_opt(j) = trace(ac.Y_k(k)*Wmus_opt*P_mmin);
    Rds_opt(j) = -trace(ac.Y_k(k)*Wmds_opt*P_mmax);
end

% normalize such that sum(d) = 1
dus_opt = normV(d_us);
dds_opt = normV(d_ds); % downspinning = upspinning in this formulation

% create scenario W_s for every mismatch
Ws_opt = zeros(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
   Ws_opt(:,:,i) = Wf_opt + Wmus_opt * min(wind.P_m(t,i), 0) + Wmds_opt * max(wind.P_m(t,i), 0);
end

decided_vars = {Wf_opt, Ws_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};

% calculate R for every scenario
R = zeros(ac.N_G, N);
% R2 = zeros(ac.N_G, N);
for i = 1:N
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R(j, i) = trace(ac.Y_k(k) * (Ws_opt(:,:,i)-Wf_opt)) - ac.C_w(k)*wind.P_m(t, i);
%         R(j, i) = trace(ac.Y_k(k) * (Wm_opt * wind.P_m(t,i))) - ac.C_w(k)*wind.P_m(t, i);
    end
end

% output ranks and table with results
format_result(ac, wind, t, decided_vars, R);

if diagnostics.problem ~= 0
    fprintf('%s (!) \t\t in %g seconds\n\n', diagnostics.info, diagnostics.solvertime);
else
    fprintf('%s \t\tin %g seconds\n\n',diagnostics.info, diagnostics.solvertime);
end

%% Check constraints for all scenarios
tol = 1e-5;
no_failed = 0;
for k = 1:ac.N_b
    no_failed = no_failed + double(ac.P_min(k)-tol > trace(ac.Y_k(k)*Wf_opt) + ac.P_D(t,k) - ac.C_w(k)*wind.P_wf(t));
    no_failed = no_failed + double(ac.P_max(k)+tol < trace(ac.Y_k(k)*Wf_opt) + ac.P_D(t,k) - ac.C_w(k)*wind.P_wf(t));
    no_failed = no_failed + double(ac.Q_min(k)-tol > trace(ac.Ybar_k(k)*Wf_opt) + ac.Q_D(t, k));
    no_failed = no_failed + double(ac.Q_max(k)+tol < trace(ac.Ybar_k(k)*Wf_opt) + ac.Q_D(t, k));
    no_failed = no_failed + double(ac.V_min(k)^2-tol > trace(ac.M_k(k)*Wf_opt));
    no_failed = no_failed + double(ac.V_max(k)^2+tol < trace(ac.M_k(k)*Wf_opt));
    no_failed = no_failed + double(abs(Wf_opt(k_ref, k_ref)) > tol);
    
    for i = 1:N
        no_failed = no_failed + double(ac.P_min(k)-tol > trace(ac.Y_k(k)*Ws_opt(:, :, i)) + ac.P_D(t,k) - ac.C_w(k)*wind.P_w(t));
        no_failed = no_failed + double(ac.P_max(k)+tol < trace(ac.Y_k(k)*Ws_opt(:, :, i)) + ac.P_D(t,k) - ac.C_w(k)*wind.P_w(t));
        no_failed = no_failed + double(ac.Q_min(k)-tol > trace(ac.Ybar_k(k)*Ws_opt(:,:,i)) + ac.Q_D(t, k));
        no_failed = no_failed + double(ac.Q_max(k)+tol < trace(ac.Ybar_k(k)*Ws_opt(:,:,i)) + ac.Q_D(t, k));
        no_failed = no_failed + double(ac.V_min(k)^2-tol > trace(ac.M_k(k)*Ws_opt(:,:,i)));
        no_failed = no_failed + double(ac.V_max(k)^2+tol < trace(ac.M_k(k)*Ws_opt(:,:,i)));
        no_failed = no_failed + double(abs(Ws_opt(k_ref, k_ref, i)) > tol);
    end
end

if ~no_failed
    fprintf('All constraints satisfied\n')
else
    fprintf('%i constraints violated\n', no_failed)
end
        