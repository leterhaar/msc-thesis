%% solve OPF-RS problem for a longer time horizon in one shot
% 24-01-17

clear
clc
yalmip('clear');
clf
ops = sdpsettings;
% ops = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1);
% ops = sdpsettings('solver', 'sedumi', 'verbose', 1, 'debug', 1);

% load models
T = 24;
N = 5;
tol = 1e-5;
ac = AC_model('case30'); 
ac.set_WPG_bus(22);

wind = wind_model(ac, 24, 0.2);
wind.generate(N);
wind.shorter_horizon(T);
% wind.plot;

%% C-OPF-RS
fprintf('RUNNING C-OPF-RS\n');

% define variables
W_f = cell(T, 1);
W_s = cell(T, N);
R_us = cell(T, 1);
R_ds = cell(T, 1);
d_us = cell(T, 1);
d_ds = cell(T, 1);

for t = 1:T
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
for t = 1:T
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
for t = 1:T
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
end

% MC simulations
wind_MC = wind_model(ac, 24, 0.2);
wind_MC.generate(1e4);
wind_MC.shorter_horizon(T);
[violations, loadings] = simulate(ac, values_cell(W_f), values_cell(d_us), values_cell(d_ds), wind_MC);

% plot
initfig('Violations', 1);
hold off
b = bar3(violations);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
ylabel('Line number');
zlabel('Violated [%]');
xlabel('Hour');
% ylim([0.5 ac.N_l+0.5]);


initfig('Loadings', 2);
boxplot(loadings');
hold on;
xlims = xlim;
plot([-1 ac.N_l+1], [100 100], '--');
xlim(xlims);
xlabel('Line number');
ylabel('Loading [%]');

%% C-OPF-RS 2
fprintf('RUNNING C-OPF-RS 2 \n');

wind2 = copy(wind);
wind2.get_extremes();

% define variables
W_f = cell(T, 1);
W_s = cell(T, 2);
R_us = cell(T, 1);
R_ds = cell(T, 1);
d_us = cell(T, 1);
d_ds = cell(T, 1);

for t = 1:T
    W_f{t} = sdpvar(2*ac.N_b);
    R_us{t} = sdpvar(ac.N_G, 1);
    R_ds{t} = sdpvar(ac.N_G, 1);
    d_us{t} = sdpvar(ac.N_G, 1);
    d_ds{t} = sdpvar(ac.N_G, 1);
    for i = 1:2
        W_s{t, i} = sdpvar(2*ac.N_b);
    end
end

Obj = 0;
C = [];

% loop over time
for t = 1:T
    % define objective
    Obj = Obj + objective(W_f{t}, R_us{t}, R_ds{t});
    
    % define constraints
    C = [C, feasibleW(W_f{t}, wind2.P_wf(t))];
    C = [C, W_f{t} >= 0];

    for i = 1:2
        C = [C, feasibleW(W_s{t, i}, wind2.P_wf(t) + wind2.P_mpos(t,i) + wind2.P_mneg(t, i))];
        C = [C, W_s{t, i} >= 0];
        
        for j = 1:ac.N_G
            k = ac.Gens(j);
            
            % Bound R between R_us and R_ds
            C = [C, -R_ds{t}(j) ...
                 <= trace(ac.Y_k(k)*(W_s{t, i}-W_f{t})) - ac.C_w(k)*(wind2.P_mpos(t, i) + wind2.P_mneg(t,i)) <= ...
                    R_us{t}(j)];
            % relate W_s and W_f through d_ds and d_us
            C = [C, trace(ac.Y_k(k)*(W_s{t, i}-W_f{t})) - ac.C_w(k)*(wind2.P_mpos(t, i) + wind2.P_mneg(t,i)) == ...
                d_us{t}(j) * wind2.P_mneg(t, i) - d_ds{t}(j) * wind2.P_mpos(t,i)];
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

% simulate_network(ac, wind2, values_cell(W_f), values_cell(d_us), values_cell(d_ds), [], 1)
%% 
% MC simulations
wind_MC = wind_model(ac, 24, 0.2);
wind_MC.generate(1e4);
wind_MC.shorter_horizon(T);
%
[violations, loadings] = simulate(ac, values_cell(W_f), values_cell(d_us), values_cell(d_ds), wind_MC, 1, 1e-2);

% plot
initfig('Violations', 1);
hold off
b = bar3(violations);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
ylabel('Line number');
zlabel('Violated [%]');
xlabel('Hour');
zlims = zlim;
zlim([0 min(100, zlims(2))])

initfig('Loadings', 2);
boxplot(loadings');
hold on;
xlims = xlim;
plot([-1 ac.N_l+1], [100 100], '--');
xlim(xlims);
xlabel('Line number');
ylabel('Loading [%]');

%% PC-OPF-RS 1
yalmip('clear');
fprintf('\nRUNNING PC-OPF-RS 1\n');

% define SDPvars
W_f = cell(T, 1);
W_us = cell(T, 1);
W_ds = cell(T, 1);

for t = 1:T
    W_f{t} = sdpvar(2*ac.N_b);
    W_us{t} = sdpvar(2*ac.N_b);
    W_ds{t} = sdpvar(2*ac.N_b);
    R_us{t} = sdpvar(ac.N_G,1);
    R_ds{t} = sdpvar(ac.N_G,1);
end

Ysum = ac.Y_sum;

Obj = 0;
C = [];
W_f = cell(T, 1);
W_s = cell(T, N);
for t = 1:T
    % define objective
    W_f{t} = W_f{t} + wind.P_wf(t)*W_us{t};
    Obj = Obj + objective(W_f{t}, R_us{t}, R_ds{t}); 
    
    % define constraints
    C = [C, feasibleW(W_f{t}, wind.P_wf)];

    C = [C, W_f{t} >= 0];
    C = [C, W_us{t} >= 0];
    C = [C, W_ds{t} >= 0];
    
    
    for i = 1:N
        W_s{t, i} = W_f{t} + min(wind.P_w(t, i), wind.P_wf(t))*W_us{t} + ...
                    max(wind.P_w(t, i) - wind.P_wf(t), 0)*W_ds{t};

        C = [C, feasibleW(W_s{t, i}, wind.P_w(t, i))];
        for j = 1:ac.N_G
            
            k = ac.Gens(j);
            
            % Bound R between R_us and R_ds
            C = [C, -R_ds{t}(j) ...
                 <= trace(ac.Y_k(k)*(W_s{t, i}-W_f{t})) - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us{t}(j)];
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

W_fopt = cell(T, 1);
d_usopt = cell(T, 1);
d_dsopt = cell(T, 1);

for t = 1:T
    % extract solution
    R_opt = nan(ac.N_G, N);
    PG_opt = nan(ac.N_G, 1);
    W_fopt{t} = value(W_f{t}) + wind.P_wf(t)*value(W_us{t});
    W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
    for i = 1:N
        W_sscen(:,:,i) = value(W_f{t}) + min(wind.P_w(t, i), wind.P_wf(t))*zero_for_nan(value(W_us{t})) + ...
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
    PGvsPW('P3', value(W_f{t}), value(W_us{t}), value(W_ds{t}));
end
% MC simulations
wind_MC = wind_model(ac, 24, 0.2);
wind_MC.generate(1e2);
wind_MC.shorter_horizon(T);
[violations, loadings] = simulate(ac, W_fopt, d_usopt, d_dsopt, wind_MC);

% plot
initfig('Violations', 1);
hold off
b = bar3(violations);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
ylabel('Line number');
zlabel('Violated [%]');
xlabel('Hour');
% ylim([0.5 ac.N_l+0.5]);


initfig('Loadings', 2);
boxplot(loadings');
hold on;
xlims = xlim;
plot([-1 ac.N_l+1], [100 100], '--');
xlim(xlims);
xlabel('Line number');
ylabel('Loading [%]');

%% PC-OPF-RS 2
yalmip('clear');
fprintf('\nRUNNING PC-OPF-RS 2\n');

W_f = cell(T,1);
Obj = 0;
C = [];

for t = 1:T
    W_f{t} = sdpvar(2*ac.N_b);
    Obj = objective_PG(W_f{t});
    C = [C, W_f{t} >= 0, feasibleW(W_f{t}, wind.P_wf(t))];
end

% optimize(C, Obj, ops);
W_fopt_previous = values_cell(W_f);


% define SDPvars
% W_ffixed = W_fopt_previous;
W_us = cell(T, 1);
W_ds = cell(T, 1);

for t = 1:T
    W_f{t} = sdpvar(2*ac.N_b);
    W_us{t} = sdpvar(2*ac.N_b);
    W_ds{t} = sdpvar(2*ac.N_b);
    R_us{t} = sdpvar(ac.N_G,1);
    R_ds{t} = sdpvar(ac.N_G,1);
end

Ysum = ac.Y_sum;

Obj = 0;
C = [];
W_s = cell(T, N);

for t = 1:T
    % define objective
    Obj = Obj + objective(W_ffixed{t}, R_us{t}, R_ds{t}); 
    
    % define constraints
    C = [C, feasibleW(W_f{t}, wind.P_wf)];

    C = [C, W_f{t} >= 0];
    C = [C, W_us{t} >= 0];
    C = [C, W_ds{t} >= 0];
    
    
    for i = 1:N
        W_s{t, i} = W_ffixed{t} + max(-wind.P_m(t, i), 0)*W_us{t} + ...
                    max(wind.P_m(t, i), 0)*W_ds{t};

        C = [C, feasibleW(W_s{t, i}, wind.P_w(t, i))];
        for j = 1:ac.N_G
            
            k = ac.Gens(j);
            
            % Bound R between R_us and R_ds
            C = [C, -R_ds{t}(j) ...
                 <= trace(ac.Y_k(k)*(W_s{t, i}-W_f{t})) - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us{t}(j)];
        end
       
    end

    % sum Wm = 1
    C = [C, (trace(Ysum * W_us{t}) == 1)];   
    C = [C, (trace(Ysum * W_ds{t}) == -1)];
    
    % Nonnegativity constraints on reserve bounds
    C = [C, R_us{t} >= 0, R_ds{t} >= 0];

end
% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
cost = value(Obj);

W_fopt = cell(T, 1);
d_usopt = cell(T, 1);
d_dsopt = cell(T, 1);
%
for t = 1:T
    % extract solution
    R_opt = nan(ac.N_G, N);
    PG_opt = nan(ac.N_G, 1);
    W_fopt{t} = W_ffixed{t};
    W_sscen = nan(2*ac.N_b, 2*ac.N_b, N);
    for i = 1:N
        W_sscen(:,:,i) = value(W_ffixed{t}) + max(-wind.P_m(t, i),0)*value(W_us{t}) + ...
                    max(wind.P_m(t, i), 0)*value(W_ds{t});    
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

    d_usopt{t} = -normV(d_usopt{t});
    d_dsopt{t} = normV(d_dsopt{t});


    format_result(ac, wind, t, {W_fopt{t}, W_sscen, nan(ac.N_G,1), nan(ac.N_G,1), d_usopt{t}, d_dsopt{t}}, R_opt);
    PGvsPW('P7', value(W_f{t}), value(W_us{t}), value(W_ds{t}));
end
% MC simulations
wind_MC = wind_model(ac, 24, 0.2);
wind_MC.generate(1e2);
wind_MC.shorter_horizon(T);
[violations, loadings] = simulate(ac, W_fopt, d_usopt, d_dsopt, wind_MC);

% plot
initfig('Violations', 1);
hold off
b = bar3(violations);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
ylabel('Line number');
zlabel('Violated [%]');
xlabel('Hour');
% ylim([0.5 ac.N_l+0.5]);


initfig('Loadings', 2);
boxplot(loadings');
hold on;
xlims = xlim;
plot([-1 ac.N_l+1], [100 100], '--');
xlim(xlims);
xlabel('Line number');
ylabel('Loading [%]');
