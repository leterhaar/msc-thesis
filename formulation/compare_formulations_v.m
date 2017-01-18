%% define models
clear
clc
if not(exist('AC_model', 'file'))
    addpath('../misc', '../networks', '../wind');
end

N = 10;
t = 5;
    
ac = AC_model('case14');
ac.set_WPG_bus(9);

wind = wind_model(ac, 24, 0.2);
wind.generate(N);
wind.plot(t);
% ops = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1);
ops = sdpsettings('solver', 'sedumi', 'verbose', 1, 'debug', 1);

%% Separate W for every scenario (P1)
yalmip('clear');

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

% [problem, info] = simulate_network(ac, wind, value(W_f), ...
%                                               value(d_us), value(d_ds), t);
% verify(not(problem), info);

solutions = struct('Name', 'P1 (separate W for every scenario)',...
                   'Cost', value(Obj),...
                   'W_f', value(W_f),...
                   'R_us', value(R_us),...
                   'R_ds', value(R_ds),...
                   'time', status.solvertime);
               
%% Using a fixed Wm with PSD on Ws
clear Obj C
yalmip('clear');

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
    W_s = W_f + W_us * max(0, -wind.P_m(t, i)) ...
              - W_ds * max(0, wind.P_m(t,i));
          
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
C = [C, (trace(Ysum * W_us) == 1)];   
C = [C, (trace(Ysum * W_ds) == 1)];

% optimize
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);

% [problem, info] = simulate_network(ac, wind, value(W_f), ...
%                                 normV(value(R_us)), normV(value(R_ds)), t);
% verify(not(problem), info);
solutions(2) = struct('Name', 'P2 (W_us, W_ds, PSD for every scenario)',...
                   'Cost', value(Obj),...
                   'W_f', value(W_f),...
                   'R_us', value(R_us),...
                   'R_ds', value(R_ds),...
                   'time', status.solvertime);
               
%% Using a fixed Wm but with other definition of W (only PSD on parts)

clear Obj C
yalmip('clear');

% define SDPvars
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);
W_0 = sdpvar(2*ac.N_b);
W_us = sdpvar(2*ac.N_b);
W_ds = sdpvar(2*ac.N_b);
d_ds = sdpvar(ac.N_G, 1);
d_us = sdpvar(ac.N_G, 1);

% define objective
W_f = W_0 + wind.P_wf(t)*W_us;
Obj = objective(W_f, R_us, R_ds);

% define constraints
C = [];
C = [C, feasibleW(W_f, wind.P_wf)];

C = [C, W_0 >= 0];
C = [C, W_us >= 0, W_ds >= 0];

for i = 1:N
    W_s = W_0 + min(wind.P_w(t, i), wind.P_wf(t))*W_us + ...
                max(wind.P_w(t, i) - wind.P_wf(t), 0)*W_ds;
          
    C = [C, feasibleW(W_s, wind.slice(i).P_w)];
    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        
%         R(k) = trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*wind.P_m(t, i);
        
        % Bound R between R_us and R_ds
%         C = [C, -R_us(j) ...
%              <= R(k) <= ...
%                 R_ds(j)];
%             
%         PGs(k) = trace(ac.Y_k(k)*W_s) + ac.P_D(t,k) - ac.C_w(k)*wind.P_w(t, i);
        PG(k) = trace(ac.Y_k(k)*W_f) + ac.P_D(t,k) - ac.C_w(k)*wind.P_wf(t);
%         R(k) = PGs(k) - PG(k);
        
        R(k) = d_us(j) * max(0, -wind.P_m(t, i)) - d_ds(j) * max(0, wind.P_m(t, i));
 
   
                

        
        C = [C, ac.P_min(k) <= PG(k)+R(k) <= ac.P_max(k)];
        C = [C, -R_ds(j) <= R(k) <= R_us(j)];
                
        % relate W_s and W_f through d_ds and d_us
%         C = [C, trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*wind.P_m(t, i) == ...
%             d_us(j) * max(0, -wind.P_m(t, i)) - d_ds(j) * max(0, wind.P_m(t, i))];
            
%             assign_cell({W_0, W_us, W_ds}, {random_symmetric_matrix(2*ac.N_b), ...
%                                        random_symmetric_matrix(2*ac.N_b), ...
%                                        random_symmetric_matrix(2*ac.N_b)});
%             test1 = value(max(0, wind.P_m(t, i))*W_ds + min(0, wind.P_m(t, i))*W_us);
%             test2 = value(W_s - W_f);
%             verify(all_close(test1, test2));

    end

%     C = [C, sum(R) + wind.P_m(t, i) == 0];
    
%     if wind.P_m(t, i) <= 0
%         C = [C, sum(R_us / wind.P_m(t, i)) == -1];
%     else
%         C = [C, sum(R_ds / wind.P_m(t, i)) == 1];
%     end
%             
end

% Nonnegativity constraints on reserve bounds
C = [C, R_us >= 0, R_ds >= 0, d_us >= 0, d_ds >= 0];

% reserve balancing constraints
C = [C, ones(1, ac.N_G)*d_us == 1, ones(1, ac.N_G)*d_ds == 1];

%     
% C = [C, sum(R_us / wind.P_m(t, 1)) == -1, ...
%         sum(R_ds / wind.P_m(t, 2)) == 1];

% sum Wm = 1
% C = [C, (trace(Ysum * W_us) == -1)];   
% C = [C, (trace(Ysum * W_ds) == -1)];

% reserve balancing constraints
% C = [C, ones(1, ac.N_G)*d_us == 1, ones(1, ac.N_G)*d_ds == 1];

% optimize
% Obj = 0;
status = optimize(C, Obj, ops);
verify(not(status.problem), status.info);
            
sum([solutions.R_us])
sum([solutions.R_ds])
 
dopt_us = value(d_us);
dopt_ds = value(d_ds);

[problem, info] = simulate_network(ac, wind, value(W_f), ...
                                dopt_us, dopt_ds, t);
verify(not(problem), info);
% solutions(3) = struct('Name', 'P3 (W_us, W_ds, PSD on W only)',...
%                    'Cost', value(Obj),...
%                    'W_f', value(W_f),...
%                    'R_us', value(R_us),...
%                    'R_ds', value(R_ds),...
%                    'time', status.solvertime);
%                
   