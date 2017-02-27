%% solve OPF-RS problem for a longer time horizon in one shot
% 24-01-17

clear
clc
yalmip('clear');
clf

% ops = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1);
% ops = sdpsettings('solver', 'sedumi', 'verbose', 1, 'debug', 1);
%% load models
N_t = 24;
N = 1000;
tol = 1e-5;
ac = AC_model('case30'); 
ac.make_model();
ac.set_WPG_bus(22);

wind = wind_model(ac, 24, 0.2);
wind.generate(N);
wind.shorter_horizon(N_t);
wind.plot;

% P1
fprintf('RUNNING P1\n');

% define variables
W_f = cell(N_t, 1);
W_s = cell(N_t, N);
R_us = cell(N_t, 1);
R_ds = cell(N_t, 1);
d_us = cell(N_t, 1);
d_ds = cell(N_t, 1);

for t = 1:N_t
    W_f{t} = sdpvar(2*ac.N_b);
    R_us{t} = sdpvar(ac.N_G, 1);
    R_ds{t} = sdpvar(ac.N_G, 1);
    d_us{t} = sdpvar(ac.N_G, 1);
    d_ds{t} = sdpvar(ac.N_G, 1);
    for i = 1:N
        W_s{t, i} = sdpvar(2*ac.N_b);
    end
end

Obj = 0;
C = [];

% loop over time
for t = 1:N_t
    % define objective
    Obj = Obj + objective(W_f{t}, R_us{t}, R_ds{t});
    
    % define constraints
    C = [C, feasibleW(W_f{t}, wind.P_wf(t))];
    C = [C, W_f{t} >= 0];

    for i = 1:N
        C = [C, feasibleW(W_s{t, i}, wind.P_w(t, i))];
        C = [C, W_s{t, i} >= 0];
        
        for j = 1:ac.N_G
            k = ac.Gens(j);
            
            % Bound R between R_us and R_ds
            C = [C, -R_ds{t}(j) ...
                 <= trace(ac.Y_k(k)*(W_s{t, i}-W_f{t})) - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us{t}(j)];
            % relate W_s and W_f through d_ds and d_us
            C = [C, trace(ac.Y_k(k)*(W_s{t, i}-W_f{t})) - ac.C_w(k)*wind.P_m(t, i) == ...
                d_us{t}(j) * max(0, -wind.P_m(t, i)) - d_ds{t}(j) * max(0, wind.P_m(t, i))];
        end

    end
    % Nonnegativity constraints on reserve bounds
    C = [C, R_us{t} >= 0, R_ds{t} >= 0];

    % reserve balancing constraints
    C = [C, ones(1, ac.N_G)*d_us{t} == 1, ones(1, ac.N_G)*d_ds{t} == 1];

end
% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

% extract solution per hour
for t = 1:N_t
    R_opt = nan(ac.N_G, N);
    PG_opt = nan(ac.N_G, 1);
    W_fopt = zero_for_nan(value(W_f{t}));
    W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
    for i = 1:N
        W_sscen(:, :, i) = value(W_s{t, i});

        for j = 1:ac.N_G
            k = ac.Gens(j);
            R_opt(j, i) = trace(ac.Y_k(k)*(W_sscen(:, :, i)-W_fopt)) - ac.C_w(k)*wind.P_m(t, i);
            if i == 1
                PG_opt(j) = trace(ac.Y_k(k)*W_fopt) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
            end
        end
    end
    R_usopt = value(R_us{t});
    R_dsopt = value(R_ds{t});
    d_usopt = value(d_us{t});
    d_dsopt = value(d_ds{t});


    format_result(ac, wind, t, {W_fopt, W_sscen, R_usopt, R_dsopt, d_usopt, d_dsopt}, R_opt);

    PGvsPW('P1', W_sscen, W_fopt);

    
end

simulate_network(ac, wind, values_cell(W_f), values_cell(d_us), values_cell(d_ds), [], 1, tol);

solutions(1) = struct('Name', 'P1', ...
                   'Cost', cost, ...
                   'solvertime', status.solvertime, ...
                   'solverinfo', status.info);

%% P3 different objective
yalmip('clear');
fprintf('\nRUNNING P3\n');


% define SDPvars
W_0 = cell(N_t, 1);
W_us = cell(N_t, 1);
W_ds = cell(N_t, 1);

for t = 1:N_t
    W_0{t} = sdpvar(2*ac.N_b);
    W_us{t} = sdpvar(2*ac.N_b);
    W_ds{t} = sdpvar(2*ac.N_b);
    R_us{t} = sdpvar(ac.N_G,1);
    R_ds{t} = sdpvar(ac.N_G,1);
end

Ysum = zeros_like(ac.Y_k(1));
Ysum_us = zeros_like(Ysum);
Ysum_ds = zeros_like(Ysum);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Ysum = Ysum + ac.Y_k(k);
    Ysum_us = Ysum_us + ac.c_us(j) * ac.Y_k(k);
    Ysum_ds = Ysum_ds + ac.c_ds(j) * ac.Y_k(k);
end

Obj = 0;
C = [];
W_f = cell(N_t, 1);
W_s = cell(N_t, N);
for t = 1:N_t
    % define objective
    W_f{t} = W_0{t} + wind.P_wf(t)*W_us{t};
    Obj = Obj + objective(W_f{t}, R_us{t}, R_ds{t}); 
    
    % define constraints
    C = [C, feasibleW(W_f{t}, wind.P_wf)];

    C = [C, W_0{t} >= 0];
    C = [C, W_us{t} >= 0];
    C = [C, W_ds{t} >= 0];
    
    
    for i = 1:N
        W_s{t, i} = W_0{t} + min(wind.P_w(t, i), wind.P_wf(t))*W_us{t} + ...
                    max(wind.P_w(t, i) - wind.P_wf(t), 0)*W_ds{t};

        C = [C, feasibleW(W_s{t, i}, wind.P_w(t, i))];
        for j = 1:ac.N_G
            
            k = ac.Gens(j);
            
            % Bound R between R_us and R_ds
            C = [C, -R_ds{t}(j) ...
                 <= trace(ac.Y_k(k)*(W_s{t, i}-W_f{t})) - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us{t}(j)];
            
            if i == 1
                C = [C, trace(ac.Y_(k) * W_us{t}) <= 0];
                C = [C, trace(ac.Y_(k) * W_ds{t}) <= 0];
            end
            
        end
       
    end

    % sum Wm = 1
    C = [C, (trace(Ysum * W_us{t}) == -1)];   
    C = [C, (trace(Ysum * W_ds{t}) == -1)];
    
    % Nonnegativity constraints on reserve bounds
    C = [C, R_us{t} >= 0, R_ds{t} >= 0];

end
% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

W_fopt = cell(N_t, 1);
d_usopt = cell(N_t, 1);
d_dsopt = cell(N_t, 1);

for t = 1:N_t
    % extract solution
    R_opt = nan(ac.N_G, N);
    PG_opt = nan(ac.N_G, 1);
    W_fopt{t} = value(W_0{t}) + wind.P_wf(t)*value(W_us{t});
    W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
    for i = 1:N
        W_sscen(:,:,i) = value(W_0{t}) + min(wind.P_w(t, i), wind.P_wf(t))*zero_for_nan(value(W_us{t})) + ...
                    max(wind.P_w(t, i) - wind.P_wf(t), 0)*zero_for_nan(value(W_ds{t}));    
        for j = 1:ac.N_G
            k = ac.Gens(j);
            R_opt(j,i) = trace(ac.Y_k(k)*(W_sscen(:,:,i)-W_fopt{t})) - ac.C_w(k)*wind.P_m(t, i);
            if i == 1
                PG_opt(j) = trace(ac.Y_k(k)*W_fopt{t}) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
            end
        end
    end

    for j = 1:ac.N_G
        k = ac.Gens(j);
        d_usopt{t}(j) = -trace(ac.Y_(k) * value(W_us{t}));
        d_dsopt{t}(j) = -trace(ac.Y_(k) * value(W_ds{t}));
    end

    d_usopt{t} = normV(d_usopt{t});
    d_dsopt{t} = normV(d_dsopt{t});


    format_result(ac, wind, t, {W_fopt{t}, W_sscen, nan(ac.N_G,1), nan(ac.N_G,1), d_usopt{t}, d_dsopt{t}}, R_opt);
    PGvsPW('P3', value(W_0{t}), value(W_us{t}), value(W_ds{t}));
end
    
simulate_network(ac, wind, W_fopt, d_usopt, d_dsopt, [], 1, tol);
% verify(not(problem), info);
solutions(2) = struct('Name', 'P3 new cost function',...
                   'Cost', cost,...
                   ...'PGCost', objective_PG(W_fopt),...
                   ...'RCost', a_posteriori_cost(PG_opt, R_opt) - objective_PG(W_fopt),...
                   ...'APCost', a_posteriori_cost(PG_opt, R_opt),...
                   ...'W_f', value(W_f),...
                   ...'R_us', [],...
                   ...'R_ds', [],...
                   'solvertime', status.solvertime,...
                   'solverinfo', status.info);
