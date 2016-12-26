% script to test whether the AC OPF with only extremes yields the same
% solution as the original fully defined problem
clc
clear
yalmip('clear');

if not(exist('AC_f', 'file'))
    addpath('../formulation');
    addpath('../wind');
    addpath('../misc');
    addpath('../experiment');
    addpath('../networks');
end

%% load models
N_t = 24;   % optimization horizon
N = 5;      % number of scenarios used for optimization
t = 1; % timestep used for this demonstration (todo: add for everything)

% load network and wind models
ac = AC_model('case_ieee30a');
ac.set_WPG_bus(10);
wind = wind_model(ac, N_t, 0.9);

% generate a number of scenarios
wind.generate(N);

% add the forecast as a scenario
wind.P_w = [wind.P_wf, wind.P_w];
wind.P_m = [zeros(N_t, 1), wind.P_m];

% optimization settings
ops = sdpsettings('solver', 'mosek', 'verbose', 1);
%% define optimization variables and formulate objective

% define Xs \in S^(2n)
Ws = cell(N+1,1);
for i = 1:N+1
    Ws{i} = sdpvar(2*ac.N_b, 2*ac.N_b, 'symmetric');
end

% define R_us, R_ds and d_us, d_ds, alpha
R_us = sdpvar(ac.N_G,1);
R_ds = sdpvar(ac.N_G,1);
d_us = sdpvar(ac.N_G,1);
d_ds = sdpvar(ac.N_G,1);
alpha = sdpvar(ac.N_G,1);


%% Define objective using epigraph notation
Obj2 = ac.c_us' * R_us + ac.c_ds' * R_ds + (ones(1, ac.N_G) * alpha);

C_epigraph = [];
for j = 1:ac.N_G
    k = ac.Gens(j);
    P_G = trace(ac.Y_k(k)*Ws{1})+ac.P_D(t,k);
    C_epigraph = [C_epigraph, ([ac.c_li(j) * P_G - alpha(j), sqrt(ac.c_qu(j)) * P_G;
                  sqrt(ac.c_qu(j)) * P_G, -1] <= 0):...
                  sprintf('Epigraph g%2i', k)];
end

%% define cost function
Obj = 0;

% real power generation in scenario 1 (= forecast)
for j = 1:ac.N_G
    k = ac.Gens(j);
    Obj = Obj + ac.c_qu(j)*((trace(Ws{1}*ac.Y_k(k))+ac.P_D(t,k))^2) ...
              + ac.c_li(j)*(trace(Ws{1}*ac.Y_k(k))+ac.P_D(t,k));
end

% add up and downspinning 
Obj = Obj + ac.c_us' * R_us + ac.c_ds' * R_ds;
%% define constraints (except psd)

C = [R_us >= 0, R_ds >= 0, sum(d_ds) == 1, sum(d_us) == 1];

% loop over scenario constraints
for i = 1:N+1
    % refbus angle constraints
    refbus_index = ac.refbus + ac.N_b;
    C = [C; Ws{i}(refbus_index, refbus_index) == 0];
            
    % psd constraints
    C = [C; Ws{i} >= 0];
    
    for k = 1:ac.N_b
        % real power injection limits
        C = [C; (ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i) <= ...
                trace(Ws{i} * ac.Y_k(k)) <= ...
                ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i)):...
                sprintf('Pinj | s%2i | b%2i', i, k)];

        % reactive power injection limits
        C = [C; (ac.Q_min(k) - ac.Q_D(t, k) <= ...
                trace(Ws{i} * ac.Ybar_k(k)) <= ...
                ac.Q_max(k) - ac.Q_D(t, k)):...
                sprintf('Qinj | s%2i | b%2i', i, k)];
        
        % voltage magnitude limits
        C = [C; ((ac.V_min(k))^2 <= ...
                trace(Ws{i} * ac.M_k(k)) <= ...
                (ac.V_max(k))^2):...
                sprintf('Vbus | s%2i | b%2i', i, k)];
    end
    
    if i > 1
        % reserve balancing constraints
        for j = 1:ac.N_G
                    % bus index
            k = ac.Gens(j);

            % Bound R between R_us and R_ds
            C = [C; (-R_ds(j) <= ...
                    trace((Ws{i} - Ws{1}) * ac.Y_k(k)) ...
                    - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us(j)):...
                    sprintf('Rdus | s%2i | b%2i', i, k)];

            % relate W_s and W_f through d_ds and d_us
            C = [C; (trace((Ws{i} - Ws{1}) * ac.Y_k(k)) ...
                    - ac.C_w(k)*wind.P_m(t, i) == ...
                    d_us(j) * max(0, -wind.P_m(t, i)) ...
                    - d_ds(j) * max(0, wind.P_m(t, i))):...
                    sprintf('Rbal | s%2i | b%2i', i, k)];        
        end
    end
    
end
% solve problem with psd constraint on every scenario matrix variable
tic
status = optimize(C, Obj, ops);
toc
verify(not(status.problem), status.info);
verify(not(any(check(C) < -1e-6)), 'Infeasible solution!');

Wsstar_all = values_cell(Ws);
Rus_all = value(R_us);
Rds_all = value(R_ds);
dus_all = value(d_us);
dds_all = value(d_ds);
alpha_all = value(alpha);
Obj_all = value(Obj2);

for i = 1:N+1
    verify(svd_rank(Wsstar_all{i}) == 1, 'X for s%i not rank 1', i);
end

%% discard the intermediate scenarios, formulate again and solve
yalmip('clear');
N = 2;
wind.use_extremes(t);
wind.P_w = [wind.P_wf, wind.P_w];
wind.P_m = [zeros(N_t, 1), wind.P_m];

% define Xs \in S^(2n)
Ws = cell(N+1,1);
for i = 1:N+1
    Ws{i} = sdpvar(2*ac.N_b, 2*ac.N_b, 'symmetric');
end

% define R_us, R_ds and d_us, d_ds, alpha
R_us = sdpvar(ac.N_G,1);
R_ds = sdpvar(ac.N_G,1);
d_us = sdpvar(ac.N_G,1);
d_ds = sdpvar(ac.N_G,1);
alpha = sdpvar(ac.N_G,1);

%% define cost function
Obj = 0;

% real power generation in scenario 1 (= forecast)
for j = 1:ac.N_G
    k = ac.Gens(j);
    Obj = Obj + ac.c_qu(j)*((trace(Ws{1}*ac.Y_k(k))+ac.P_D(t,k))^2) ...
              + ac.c_li(j)*(trace(Ws{1}*ac.Y_k(k))+ac.P_D(t,k));
end

% add up and downspinning 
Obj = Obj + ac.c_us' * R_us + ac.c_ds' * R_ds;

%% define constraints (for extremes only)

C = [R_us >= 0, R_ds >= 0, sum(d_ds) == 1, sum(d_us) == 1];

% loop over scenario constraints
for i = 1:N+1
    % refbus angle constraints
    refbus_index = ac.refbus + ac.N_b;
    C = [C; Ws{i}(refbus_index, refbus_index) == 0];
            
    % psd constraints
    C = [C; Ws{i} >= 0];
    
    for k = 1:ac.N_b
        % real power injection limits
        C = [C; (ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i) <= ...
                trace(Ws{i} * ac.Y_k(k)) <= ...
                ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i)):...
                sprintf('Pinj | s%2i | b%2i', i, k)];

        % reactive power injection limits
        C = [C; (ac.Q_min(k) - ac.Q_D(t, k) <= ...
                trace(Ws{i} * ac.Ybar_k(k)) <= ...
                ac.Q_max(k) - ac.Q_D(t, k)):...
                sprintf('Qinj | s%2i | b%2i', i, k)];
        
        % voltage magnitude limits
        C = [C; ((ac.V_min(k))^2 <= ...
                trace(Ws{i} * ac.M_k(k)) <= ...
                (ac.V_max(k))^2):...
                sprintf('Vbus | s%2i | b%2i', i, k)];
    end
    
    if i > 1
        % reserve balancing constraints
        for j = 1:ac.N_G
                    % bus index
            k = ac.Gens(j);

            % Bound R between R_us and R_ds
            C = [C; (-R_ds(j) <= ...
                    trace((Ws{i} - Ws{1}) * ac.Y_k(k)) ...
                    - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us(j)):...
                    sprintf('Rdus | s%2i | b%2i', i, k)];

            % relate W_s and W_f through d_ds and d_us
            C = [C; (trace((Ws{i} - Ws{1}) * ac.Y_k(k)) ...
                    - ac.C_w(k)*wind.P_m(t, i) == ...
                    d_us(j) * max(0, -wind.P_m(t, i)) ...
                    - d_ds(j) * max(0, wind.P_m(t, i))):...
                    sprintf('Rbal | s%2i | b%2i', i, k)];        
        end
    end
    
end
%% solve problem with psd constraint on every scenario matrix variable
tic
status = optimize(C, Obj, ops);
toc
verify(not(status.problem), status.info);
verify(not(any(check(C) < -1e-6)), 'Infeasible solution!');

Wsstar_extremes = values_cell(Ws);
Rus_extremes = value(R_us);
Rds_extremes = value(R_ds);
dus_extremes = value(d_us);
dds_extremes = value(d_ds);
alpha_extremes = value(alpha);
Obj_extremes = value(Obj);

for i = 1:N+1
    verify(svd_rank(Wsstar_extremes{i}) == 1, 'X for s%i not rank 1', i);
    verify(is_psd(Wsstar_extremes{i}), 'W not psd for s%2i', i);
end

verify(all_close(Wsstar_extremes{1}, Wsstar_all{1}), 'Wf not the same: %4e', norm(Wsstar_extremes{1}-Wsstar_all{1}, 'inf'));

% extract and compare voltage vectors 
X_all = sqrt(diag(Wsstar_all{1})) .* sign(Wsstar_all{1}(:,1));
X_extremes = sqrt(diag(Wsstar_extremes{1})) .* sign(Wsstar_extremes{1}(:,1));
verify(all_close(X_all, X_extremes), 'Voltage vectors not close: %4e', norm(X_all-X_extremes, 'inf'));

% extract PG for all busses 
Pg_all = zeros(ac.N_b,1);
Pg_extremes = zeros_like(Pg_all);

for k = 1:ac.N_b
    Pg_all(k) = trace(ac.Y_k(k)*Wsstar_all{1}) + ac.P_D(t, k) - ac.C_w(k) * wind.P_wf(t);
    Pg_extremes(k) = trace(ac.Y_k(k)*Wsstar_extremes{1}) + ac.P_D(t, k) - ac.C_w(k) * wind.P_wf(t);
end

% verify(all_close(Pg_all, Pg_extremes), 'Pg not close: %4e', norm(Pg_all - Pg_extremes, 'inf'));
% verify(all_close(Rus_extremes, Rus_all, 1e-5), 'Rus not close: %4e', norm(Rus_extremes-Rus_all, 'inf'));
% verify(all_close(Rds_extremes, Rds_all, 1e-5), 'Rds not close: %4e', norm(Rds_extremes-Rds_all, 'inf'));
% verify(all_close(dus_extremes, dus_all, 1e-4), 'dus not close: %4e', norm(dus_extremes-dus_all, 'inf'));
% verify(all_close(dds_extremes, dds_all, 1e-4), 'dds not close: %4e', norm(dds_extremes-dds_all, 'inf'));

Wf_all = Wsstar_all{1};
Wf_extremes = Wsstar_extremes{1};
%%
var_names = {'Wf', 'Pg', 'X', 'Rus', 'Rds', 'dus', 'dds'};
i = 1;

for i = 1:length(var_names);   
    var_name = var_names{i};
    fprintf('%s & %e\n',var_name, ...
       eval(['norm(' var_name '_all - ' var_name '_extremes, ''inf'');']));

end
