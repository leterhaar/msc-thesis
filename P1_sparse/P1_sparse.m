%% P1 sparse
% runs the formulation P1 with the PSD constraints on submatrices based on
% the tree decomposition


%% initiate models
clear
yalmip('clear');

addpath('../experiments');
addpath('../formulation');
addpath('../misc');

exp = experiment.load('debug14-100')

N = exp.wind.N;
T = exp.wind.horizon;
ac = exp.network.model;
ac.get_bags();
wind = exp.wind.model;

%% define problem
W_f = cell(T, 1);
W_s = cell(T, N);
R_us = cell(T, 1);
R_ds = cell(T, 1);
d_us = cell(T, 1);
d_ds = cell(T, 1);

C = []; Obj = 0;

p = progress('Definining problem', (N+1)*T);
for t = 1:T

    % define SDPvars
    W_f{t} = sdpvar(2*ac.N_b);
    for i = 1:N
        W_s{t, i} = sdpvar(2*ac.N_b);
    end

    R_us{t} = sdpvar(ac.N_G, 1);
    R_ds{t} = sdpvar(ac.N_G, 1);
    d_us{t} = sdpvar(ac.N_G, 1);
    d_ds{t} = sdpvar(ac.N_G, 1);

    % define objective 
    Obj = Obj + objective(W_f{t}, R_us{t}, R_ds{t});
    
    % feasible and PSD Wf
    C = [C, feasibleW(W_f{t}, wind.P_wf(t))];
    C = [C, PSD_bags(W_f{t}, ac.bags)];
    
    % nonnegativity Rus, Rds
    C = [C, R_us{t} >= 0, R_ds{t} >= 0];
    
    % sum d == 1
    C = [C, sum(d_us{t}) == 1, sum(d_ds{t}) == 1];
    p.ping();

    for i = 1:N
        
        % feasible and PSD Ws
        C = [C, feasibleW(W_s{t, i}, wind.P_w(t, i))];
        C = [C, PSD_bags(W_s{t, i}, ac.bags)];
        
        % bound reserve
        for j = 1:ac.N_G
            k = ac.Gens(j);
            C = [C, -R_ds{t}(j) <= ...
                    trace(ac.Y_(k) * (W_s{t, i} - W_f{t})) - ac.C_w(k)*wind.P_m(t, i) ...
                    <= R_us{t}(j)];
            C = [C, trace(ac.Y_(k) * (W_s{t, i} - W_f{t})) - ac.C_w(k)*wind.P_m(t, i) ...
                == -d_us{t}(j) * min(0, wind.P_m(t, i)) - d_ds{t}(j) * max(0, wind.P_m(t, i))]; 
        end
        p.ping()
        
    end
end

%% solve problem
status = optimize(C, Obj)

residuals = check(C);

%% store solutions
R_usopt = values_cell(R_us);
R_dsopt = values_cell(R_ds);
d_usopt = values_cell(d_us);
d_dsopt = values_cell(d_ds);
W_fopt = values_cell(W_f);
W_sopt = values_cell(W_s);

