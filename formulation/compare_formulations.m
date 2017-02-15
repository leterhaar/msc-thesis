%% define models
clear
clc
if not(exist('AC_model', 'file'))
    addpath('../misc', '../networks', '../wind');
end

N = 2;
t = 15;
tol = 1e-4;


ac = AC_model('case14');
ac.set_WPG_bus(9);
% ac.c_ds = 1./ac.c_us;

wind = wind_model(ac, 24, 0.2);
while 1
    wind.generate(N);
    wind.plot(t);
    pause(0.01);
    if min(wind.P_w(t,:)) == 0
        break
    end
end
% wind.use_extremes(t);
% N = 2;
wind.plot(t);
ops = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1);
% ops = sdpsettings('solver', 'sedumi', 'verbose', 1, 'debug', 1);

%% P1
yalmip('clear');
fprintf('\nRUNNING P1\n');

% define SDPvars
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);
d_us = sdpvar(ac.N_G, 1);
d_ds = sdpvar(ac.N_G, 1);
W_f = sdpvar(2*ac.N_b);
W_s = cell(N,1);

% define objective
Obj = objective(W_f, R_us, R_ds);

% define constraints
C = feasibleW(W_f, wind.P_wf);
C = [C, W_f >= 0];

for i = 1:N
    W_s{i} = sdpvar(2*ac.N_b);
    C = [C, feasibleW(W_s{i}, wind.slice(i).P_w)];
    C = [C, W_s{i} >= 0];
    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        
        % Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= trace(ac.Y_k(k)*(W_s{i}-W_f)) - ac.C_w(k)*wind.P_m(t, i) <= ...
                R_us(j)];
        % relate W_s and W_f through d_ds and d_us
        C = [C, trace(ac.Y_k(k)*(W_s{i}-W_f)) - ac.C_w(k)*wind.P_m(t, i) == ...
            d_us(j) * max(0, -wind.P_m(t, i)) - d_ds(j) * max(0, wind.P_m(t, i))];
    end

end

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];

% reserve balancing constraints
C = [C, ones(1, ac.N_G)*d_us == 1, ones(1, ac.N_G)*d_ds == 1];

% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

R_opt = nan(ac.N_G, N);
PG_opt = nan(ac.N_G, 1);
W_fopt = zero_for_nan(value(W_f));
W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    W_sscen(:,:,i) = value(W_s{i});
          
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R_opt(j,i) = trace(ac.Y_k(k)*(W_sscen(:,:,i)-W_fopt)) - ac.C_w(k)*wind.P_m(t, i);
        if i == 1
            PG_opt(j) = trace(ac.Y_k(k)*W_fopt) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
        end
    end
end
R_usopt = value(R_us);
R_dsopt = value(R_ds);
d_usopt = value(d_us);
d_dsopt = value(d_ds);


format_result(ac, wind, t, {W_fopt, W_sscen, R_usopt, R_dsopt, d_usopt, d_dsopt}, R_opt);

PGvsPW('P1', W_sscen, W_fopt);

simulate_network(ac, wind, value(W_f), value(d_us), value(d_ds), t, 1, tol);


solutions(1) = struct('Name', 'P1',...
                   'Cost', cost,...
                   'PGCost', objective_PG(W_fopt),...
                   'RCost', a_posteriori_cost(PG_opt, R_opt) - objective_PG(W_fopt),...
                   'APCost', a_posteriori_cost(PG_opt, R_opt),...
                   'W_f', value(W_f),...
                   'R_us', value(R_us),...
                   'R_ds', value(R_ds),...
                   'solvertime', status.solvertime,...
                   'solverinfo', status.info);
               
%% P1 different objective
yalmip('clear');
fprintf('\nRUNNING P1 different objective\n');

% define SDPvars
d_us = sdpvar(ac.N_G, 1);
d_ds = sdpvar(ac.N_G, 1);
W_f = sdpvar(2*ac.N_b);
W_s = cell(N,1);

Ysum = zeros_like(ac.Y_k(1));
Ysum_us = zeros_like(Ysum);
Ysum_ds = zeros_like(Ysum);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Ysum = Ysum + ac.Y_k(k);
    Ysum_us = Ysum_us + ac.c_us(j) * ac.Y_k(k);
    Ysum_ds = Ysum_ds + ac.c_ds(j) * ac.Y_k(k);
end


% define objective
Obj = objective_PG(W_f);

% define constraints
C = feasibleW(W_f, wind.P_wf);
C = [C, W_f >= 0];

for i = 1:N
    W_s{i} = sdpvar(2*ac.N_b);
    C = [C, feasibleW(W_s{i}, wind.slice(i).P_w)];
    C = [C, W_s{i} >= 0];
    
   
    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        
        % relate W_s and W_f through d_ds and d_us
        C = [C, trace(ac.Y_k(k)*(W_s{i}-W_f)) - ac.C_w(k)*wind.P_m(t, i) == ...
            d_us(j) * max(0, -wind.P_m(t, i)) - d_ds(j) * max(0, wind.P_m(t, i))];
        
        % add to objective
        Obj = Obj + (2/N)*(ac.c_us(j) * d_us(j) * max(0, -wind.P_m(t, i)) + ... 
              ac.c_ds(j) * d_ds(j) * max(0, wind.P_m(t, i)));
    end

end

% reserve balancing constraints
C = [C, ones(1, ac.N_G)*d_us == 1, ones(1, ac.N_G)*d_ds == 1];
C = [C, d_us >= 0, d_ds >= 0];

% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

R_opt = nan(ac.N_G, N);
PG_opt = nan(ac.N_G, 1);
W_fopt = zero_for_nan(value(W_f));
W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    W_sscen(:,:,i) = value(W_s{i});
          
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R_opt(j,i) = trace(ac.Y_k(k)*(W_sscen(:,:,i)-W_fopt)) - ac.C_w(k)*wind.P_m(t, i);
        if i == 1
            PG_opt(j) = trace(ac.Y_k(k)*W_fopt) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
        end
    end
end
R_usopt = zeros(ac.N_G, 1);
R_dsopt = zeros(ac.N_G, 1);
d_usopt = value(d_us);
d_dsopt = value(d_ds);


format_result(ac, wind, t, {W_fopt, W_sscen, R_usopt, R_dsopt, d_usopt, d_dsopt}, R_opt);

PGvsPW('P1', W_sscen, W_fopt);

[problem, info] = simulate_network(ac, wind, value(W_f), ...
                                              value(d_us), value(d_ds), t, 1, tol);
% verify(not(problem), info);

solutions(2) = struct('Name', 'P1 new cost function',...
                   'Cost', cost,...
                   'PGCost', objective_PG(W_fopt),...
                   'RCost', a_posteriori_cost(PG_opt, R_opt) - objective_PG(W_fopt),...
                   'APCost', a_posteriori_cost(PG_opt, R_opt),...
                   'W_f', value(W_f),...
                   'R_us', [],...
                   'R_ds', [],...
                   'solvertime', status.solvertime,...
                   'solverinfo', status.info);
               

%% P2

clear Obj C
yalmip('clear');
fprintf('\nRUNNING P2\n');


% define SDPvars
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);
W_f = sdpvar(2*ac.N_b);
W_us = sdpvar(2*ac.N_b);
W_ds = sdpvar(2*ac.N_b);

% define objective
Obj = objective(W_f, R_us, R_ds);

% define constraints
C = feasibleW(W_f, wind.P_wf);
C = [C, W_f >= 0];

for i = 1:N
    W_s = W_f + W_us * max(0, wind.P_m(t, i)) ...
              + W_ds * min(0, wind.P_m(t,i));
          
    C = [C, feasibleW(W_s, wind.slice(i).P_w)];
    C = [C, W_s >= 0];
    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        
        % Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*wind.P_m(t, i) <= ...
                R_us(j)];
    end

end

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];

% reserve balancing constraints
Ysum = zeros_like(ac.Y_k(1));
for k = ac.Gens'
    Ysum = Ysum + ac.Y_k(k);
end

% sum Wm = 1
C = [C, (trace(Ysum * W_us) == -1)];   
C = [C, (trace(Ysum * W_ds) == -1)];

% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

R_opt = nan(ac.N_G, N);
PG_opt = nan(ac.N_G, 1);
W_fopt = zero_for_nan(value(W_f));
W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    W_sscen(:,:,i) =  W_fopt + value(W_us) * max(0, wind.P_m(t, i)) ...
              + value(W_ds) * min(0, wind.P_m(t,i));
          
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R_opt(j,i) = trace(ac.Y_k(k)*(W_sscen(:,:,i)-W_fopt)) - ac.C_w(k)*wind.P_m(t, i);
        if i == 1
            PG_opt(j) = trace(ac.Y_k(k)*W_fopt) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
        end
    end
end
R_usopt = value(R_us);
R_dsopt = value(R_ds);


for j = 1:ac.N_G
    k = ac.Gens(j);
    d_usopt(j) = -trace(ac.Y_(k) * value(W_us));
    d_dsopt(j) = -trace(ac.Y_(k) * value(W_ds));
end
d_usopt = normV(d_usopt);
d_dsopt = normV(d_dsopt);


format_result(ac, wind, t, {W_fopt, W_sscen, R_usopt, R_dsopt, d_usopt, d_dsopt}, R_opt);


[problem, info] = simulate_network(ac, wind, W_fopt, d_usopt, d_dsopt, t, 1, tol);
% verify(not(problem), info);
solutions(3) = struct('Name', 'P2',...
                   'Cost', cost,...
                   'PGCost', objective_PG(W_fopt),...
                   'RCost', a_posteriori_cost(PG_opt, R_opt) - objective_PG(W_fopt),...
                   'APCost', a_posteriori_cost(PG_opt, R_opt),...
                   'W_f', value(W_f),...
                   'R_us', value(R_us),...
                   'R_ds', value(R_ds),...
                   'solvertime', status.solvertime,...
                   'solverinfo', status.info);

%% P3

clear Obj C
yalmip('clear');
fprintf('\nRUNNING P3\n');


% define SDPvars
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);
W_0 = sdpvar(2*ac.N_b);
W_us = sdpvar(2*ac.N_b);
W_ds = sdpvar(2*ac.N_b);

Ysum = zeros_like(ac.Y_k(1));
for j = 1:ac.N_G
    k = ac.Gens(j);
    Ysum = Ysum + ac.Y_k(k);
end


% define objective
W_f = W_0 + wind.P_wf(t)*W_us;
Obj = objective(W_f, R_us, R_ds);

% define constraints
C = [];
C = [C, feasibleW(W_f, wind.P_wf)];

C = [C, W_0 >= 0];
C = [C, W_us >= 0];
C = [C, W_ds >= 0];
W_s = cell(N,1);
for i = 1:N
    W_s{i} = W_0 + min(wind.P_w(t, i), wind.P_wf(t))*W_us + ...
                max(wind.P_w(t, i) - wind.P_wf(t), 0)*W_ds;

    C = [C, feasibleW(W_s{i}, wind.slice(i).P_w)];
    
    for j = 1:ac.N_G
        
        k = ac.Gens(j);
        
        % Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= trace(ac.Y_k(k)*(W_s{i}-W_f)) - ac.C_w(k)*wind.P_m(t, i) <= ...
                R_us(j)];  
        if i == 1
            C = [C, trace(ac.Y_(k) * W_us) <= 0];
            C = [C, trace(ac.Y_(k) * W_ds) <= 0];
        end
        
    end
    
end

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];

% sum Wm = 1
C = [C, (trace(Ysum * W_us) == -1)];   
C = [C, (trace(Ysum * W_ds) == -1)];

% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

% extract solution
R_opt = nan(ac.N_G, N);
PG_opt = nan(ac.N_G, 1);
W_fopt = zero_for_nan(value(W_f));
W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    W_sscen(:,:,i) = value(W_0) + min(wind.P_w(t, i), wind.P_wf(t))*zero_for_nan(value(W_us)) + ...
                max(wind.P_w(t, i) - wind.P_wf(t), 0)*zero_for_nan(value(W_ds));    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R_opt(j,i) = trace(ac.Y_k(k)*(W_sscen(:,:,i)-W_fopt)) - ac.C_w(k)*wind.P_m(t, i);
        if i == 1
            PG_opt(j) = trace(ac.Y_k(k)*W_fopt) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
        end
    end
end
R_usopt = value(R_us);
R_dsopt = value(R_ds);


for j = 1:ac.N_G
    k = ac.Gens(j);
    d_usopt(j) = -trace(ac.Y_(k) * value(W_us));
    d_dsopt(j) = -trace(ac.Y_(k) * value(W_ds));
end
d_usopt = normV(d_usopt);
d_dsopt = normV(d_dsopt);


format_result(ac, wind, t, {W_fopt, W_sscen, R_usopt, R_dsopt, d_usopt, d_dsopt}, R_opt);
PGvsPW('P3', value(W_0), value(W_us), value(W_ds));

[problem, info] = simulate_network(ac, wind, W_fopt, d_usopt, d_dsopt, t, 1, tol);
% verify(not(problem), info);
solutions(4) = struct('Name', 'P3',...
                   'Cost', cost,...
                   'PGCost', objective_PG(W_fopt),...
                   'RCost', a_posteriori_cost(PG_opt, R_opt) - objective_PG(W_fopt),...
                   'APCost', a_posteriori_cost(PG_opt, R_opt),...
                   'W_f', value(W_f),...
                   'R_us', value(R_us),...
                   'R_ds', value(R_ds),...
                   'solvertime', status.solvertime,...
                   'solverinfo', status.info);
%% P3 different objective

clear Obj C
yalmip('clear');
fprintf('\nRUNNING P3 with new cost function\n');


% define SDPvars
W_0 = sdpvar(2*ac.N_b);
W_us = sdpvar(2*ac.N_b);
W_ds = sdpvar(2*ac.N_b);

Ysum = zeros_like(ac.Y_k(1));
Ysum_us = zeros_like(Ysum);
Ysum_ds = zeros_like(Ysum);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Ysum = Ysum + ac.Y_k(k);
    Ysum_us = Ysum_us + ac.c_us(j) * ac.Y_k(k);
    Ysum_ds = Ysum_ds + ac.c_ds(j) * ac.Y_k(k);
end


% define objective
W_f = W_0 + wind.P_wf(t)*W_us;
Obj = 0;
Obj = Obj + objective_PG(W_f); % only for P_G

% define constraints
C = [];
C = [C, feasibleW(W_f, wind.P_wf)];

C = [C, W_0 >= 0];
C = [C, W_us >= 0];
C = [C, W_ds >= 0];

for i = 1:N
    W_s = W_0 + min(wind.P_w(t, i), wind.P_wf(t))*W_us + ...
                max(wind.P_w(t, i) - wind.P_wf(t), 0)*W_ds;

    C = [C, feasibleW(W_s, wind.slice(i).P_w)];
    for j = 1:ac.N_G
        
        k = ac.Gens(j);
        
        if i == 1
            C = [C, trace(ac.Y_(k) * W_us) <= 0];
            C = [C, trace(ac.Y_(k) * W_ds) <= 0];
        end
        
    end
    
    Obj = Obj + (2/N)*trace(Ysum_us * (min(0, wind.P_m(t,i   ))*W_us)) + ...
                (2/N)*trace(Ysum_ds * (max(0, wind.P_m(t,i))*W_ds));
    
end

% Nonnegativity constraints on reserve bounds
% C = [C, Q_us >= 0, Q_ds >= 0];

% sum Wm = 1
C = [C, (trace(Ysum * W_us) == -1)];   
C = [C, (trace(Ysum * W_ds) == -1)];

% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

% extract solution
R_opt = nan(ac.N_G, N);
PG_opt = nan(ac.N_G, 1);
W_fopt = zero_for_nan(value(W_f));
W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    W_sscen(:,:,i) = value(W_0) + min(wind.P_w(t, i), wind.P_wf(t))*zero_for_nan(value(W_us)) + ...
                max(wind.P_w(t, i) - wind.P_wf(t), 0)*zero_for_nan(value(W_ds));    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R_opt(j,i) = trace(ac.Y_k(k)*(W_sscen(:,:,i)-W_fopt)) - ac.C_w(k)*wind.P_m(t, i);
        if i == 1
            PG_opt(j) = trace(ac.Y_k(k)*W_fopt) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
        end
    end
end


for j = 1:ac.N_G
    k = ac.Gens(j);
    d_usopt(j) = -trace(ac.Y_(k) * value(W_us));
    d_dsopt(j) = -trace(ac.Y_(k) * value(W_ds));
end

d_usopt = normV(d_usopt);
d_dsopt = normV(d_dsopt);


format_result(ac, wind, t, {W_fopt, W_sscen, nan(ac.N_G,1), nan(ac.N_G,1), d_usopt, d_dsopt}, R_opt);
PGvsPW('P3', value(W_0), value(W_us), value(W_ds));

[problem, info] = simulate_network(ac, wind, W_fopt, d_usopt, d_dsopt, t, 1, tol);
% verify(not(problem), info);
solutions(5) = struct('Name', 'P3 new cost function',...
                   'Cost', cost,...
                   'PGCost', objective_PG(W_fopt),...
                   'RCost', a_posteriori_cost(PG_opt, R_opt) - objective_PG(W_fopt),...
                   'APCost', a_posteriori_cost(PG_opt, R_opt),...
                   'W_f', value(W_f),...
                   'R_us', [],...
                   'R_ds', [],...
                   'solvertime', status.solvertime,...
                   'solverinfo', status.info);
%% P4
Pwmax_shifts = 0.1:0.1:0.5;
Pwmax_shifts = 0.1;
for Pwmax_shift = Pwmax_shifts
    
clear Obj C
yalmip('clear');
fprintf('\nRUNNING P4\n');

Pwmax = max(wind.P_w(t,:))+Pwmax_shift; % maximum output based on WPG penetration

% define SDPvars
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);
W_max = sdpvar(2*ac.N_b);
W_us = sdpvar(2*ac.N_b);
W_ds = sdpvar(2*ac.N_b);

Ysum = zeros_like(ac.Y_k(1));
Ysum_us = zeros_like(Ysum);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Ysum = Ysum + ac.Y_k(k);
    Ysum_us = Ysum_us + ac.c_us(j) * ac.Y_k(k);
end

% define objective
W_f = W_max + (Pwmax - wind.P_wf(t))*W_ds;
Obj = objective(W_f, R_us, R_ds);

% define constraints
C = [];
C = [C, feasibleW(W_f, wind.P_wf)];

C = [C, W_max >= 0];
C = [C, W_us >= 0];
C = [C, W_ds >= 0];

for i = 1:N
    W_s{i} = W_max + max((-wind.P_w(t, i) + wind.P_wf(t)), 0) * W_us + ...
                  min(-wind.P_w(t, i) + Pwmax, -wind.P_wf(t) + Pwmax) * W_ds;
              
    C = [C, feasibleW(W_s{i}, wind.slice(i).P_w)];
    
    for j = 1:ac.N_G
        
        k = ac.Gens(j);

%         Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= trace(ac.Y_k(k)*(W_s{i}-W_f)) - ac.C_w(k)*wind.P_m(t, i) <= ...
                R_us(j)];  
        if i == 1
            C = [C, trace(ac.Y_(k) * W_us) >= 0];
            C = [C, trace(ac.Y_(k) * W_ds) >= 0];
        end
    end
    
end

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0];

% reserve balancing constraints

C = [C, (trace(Ysum * W_us) == 1)];   
C = [C, (trace(Ysum * W_ds) == 1)];

% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

R_opt = nan(ac.N_G, N);
PG_opt = nan(ac.N_G, 1);
W_fopt = zero_for_nan(value(W_f));
W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    W_sscen(:,:,i) =  value(W_max) + max((-wind.P_w(t, i) + wind.P_wf(t)), 0) * value(W_us) + ...
                  min(-wind.P_w(t, i) + Pwmax, -wind.P_wf(t) + Pwmax) * value(W_ds);
    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R_opt(j,i) = trace(ac.Y_k(k)*(W_sscen(:,:,i)-W_fopt)) - ac.C_w(k)*wind.P_m(t, i);
        if i == 1
            PG_opt(j) = trace(ac.Y_k(k)*W_fopt) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
        end
    end
end
R_usopt = value(R_us);
R_dsopt = value(R_ds);


for j = 1:ac.N_G
    k = ac.Gens(j);
    d_usopt(j) = trace(ac.Y_(k) * value(W_us));
    d_dsopt(j) = trace(ac.Y_(k) * value(W_ds));
end
d_usopt = normV(d_usopt);
d_dsopt = normV(d_dsopt);


format_result(ac, wind, t, {W_fopt, W_sscen, R_usopt, R_dsopt, d_usopt, d_dsopt}, R_opt);
PGvsPW('P4', value(W_max), value(W_us), value(W_ds), Pwmax);

[problem, info] = simulate_network(ac, wind, W_fopt, d_usopt, d_dsopt, t, 1, tol);
% verify(not(problem), info);
if any(check(C) < -tol)
    warning('Not feasible!');
end
solutions(6) = struct('Name', sprintf('P4 with shift %g', Pwmax_shift),...
                   'Cost', cost,...
                   'PGCost', objective_PG(W_fopt),...
                   'RCost', a_posteriori_cost(PG_opt, R_opt) - objective_PG(W_fopt),...
                   'APCost', a_posteriori_cost(PG_opt, R_opt),...
                   'W_f', value(W_f),...
                   'R_us', value(R_us),...
                   'R_ds', value(R_ds),...
                   'solvertime', status.solvertime,...
                   'solverinfo', status.info);
end