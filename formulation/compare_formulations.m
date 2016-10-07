%% COMPARE FORMULATIONS
% Script to run optimizations of different sizes network, number of
% scenarios, values of lambda for the three different networks
addpath('../misc');
addpath('../wind');
%% Define arrays to loop over
Ns = [100 250 500 1000 2000];
networks = {'case_ieee30'};
wind_bus = [9,              9,           22]; % corresponding wind bus
lambdas = 0.1:0.2:1;
t = 1;
r = 0;

results = struct(   'formulation',{},...
                    'network',{},...
                    'scenarios',{},...
                    'lambda',{},...
                    'prep_time',{},...
                    'opt_time',{},...
                    'failed',{},...
                    'mean_abs_diff',{},...
                    'rank_perc',{},...
                    'exit_code',{},...
                    'decided_vars',{});
results(length(networks)*length(Ns)).failed = 0; % preallocate results
%% TEST P1
p = formulate('P1');
lambda = 0.5; % use fixed lambda for this
% Loop over different networks
for n = 1:length(networks)
    
    % load model
    ac = AC_model(networks{n});
    ac.set_WPG_bus(wind_bus(n));
    
    
    % define wind model with 24 timesteps and 20 percent WPG penetration
    wind = wind_model(ac, 24, 0.2);
    
    % define decision variables (except W_s)
    W_f = sdpvar(2*ac.N_b); 
    R_us = sdpvar(ac.N_G, 1);
    R_ds = sdpvar(ac.N_G, 1);
    d_us = sdpvar(ac.N_G, 1);
    d_ds = sdpvar(ac.N_G, 1);
    
    % loop over different scenario sizes
    for N = Ns
        
        fprintf('==================================================\nPreparing P1 %s with %i scenarios ...', networks{n}, N);

        % define W_s
        W_s = sdpvar(2*ac.N_b, 2*ac.N_b, N);

        % generate dummy scenarios
        wind.dummy(N);

        decision_vars = {W_f, W_s, R_us, R_ds, d_us, d_ds};

        % prepare problem
        prep_start = tic;
        [Obj, Cons] = p.prepare(ac, wind, t, decision_vars, lambda);
        prep_time = toc(prep_start);
        fprintf('\b\b\b\b: done in %g seconds\r', prep_time);

        
        % optimize
        opts = sdpsettings('verbose', 0, 'solver', 'mosek');
        fprintf('Optimizing ...');
        opt_start = tic;
        diagnostics = optimize(Cons, Obj, opts);
        opt_time = toc(opt_start);
        fprintf('\b\b\b\b: done in %g seconds\r', opt_time);
        
        % a posteriori checking
        Wm_opt = value(W_s);
        Wf_opt = value(W_f);
        Rus_opt = value(R_us);
        Rds_opt = value(R_ds);
        dus_opt = normV(value(d_us));
        dds_opt = normV(value(d_ds));

        decided_vars = {Wf_opt, Wm_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};
        
        % check constraint satisfaction
        failed = p.evaluate(ac, wind, t, decided_vars);
        failed_perc = failed / length(p.checked_constraints);
        
        % calculate R and rank W_s for every scenario 
        R = zeros(ac.N_G, N);
        ranks = zeros(N+1, 1);
        difference = zeros(N, 1);
        for i = 1:N
            for j = 1:ac.N_G
                k = ac.Gens(j);
                R(j, i) = trace(ac.Y_k(k) * (Wm_opt(:,:,i)-Wf_opt)) - ac.C_w(k)*wind.P_m(t, i);
            end
            ranks(i) = svd_rank(Wm_opt(:,:,i));
            difference(i) = sum(R(:,i))+wind.P_m(t, i);
        end
        ranks(end) = svd_rank(Wf_opt);
        
        % calculate percentage of system matrices with rank 1
        rank_perc = sum(ranks==1)/(N+1)*100;
        
        % store results
        r = r + 1;
        results(r).formulation = 'P1';
        results(r).network = networks{n};
        results(r).scenarios = N;
        results(r).lambda = lambda;
        results(r).prep_time = prep_time;
        results(r).opt_time = opt_time;
        results(r).failed = failed_perc;
        results(r).mean_abs_diff = meanabs(difference);
        results(r).rank_perc = rank_perc;
        results(r).exit_code = diagnostics.problem;
        results(r).decided_vars = decided_vars;
    end
end
%% TEST P3
p = formulate('P3');
lambda = 0.5; % use fixed lambda for this

% Loop over different networks
for n = 1:length(networks)
    
    % load model
    ac = AC_model(networks{n});
    ac.set_WPG_bus(wind_bus(n));
    
    
    % define wind model with 24 timesteps and 20 percent WPG penetration
    wind = wind_model(ac, 24, 0.2);
    
    % define decision variables
    W_f = sdpvar(2*ac.N_b); 
    W_m = sdpvar(2*ac.N_b);
    R_us = sdpvar(ac.N_G, 1);
    R_ds = sdpvar(ac.N_G, 1);
    
    % loop over different scenario sizes
    for N = Ns
        
        fprintf('==================================================\nPreparing P3 %s with %i scenarios ...', networks{n}, N);

        % generate dummy scenarios
        wind.dummy(N);

        decision_vars = {W_f, W_m, R_us, R_ds};

        % prepare problem
        prep_start = tic;
        [Obj, Cons] = p.prepare(ac, wind, t, decision_vars, lambda);
        prep_time = toc(prep_start);
        fprintf('\b\b\b\b: done in %g seconds\n', prep_time);
        
        % optimize
        opts = sdpsettings('verbose', 0, 'solver', 'mosek');
        fprintf('Optimizing ...');
        opt_start = tic;
        diagnostics = optimize(Cons, Obj, opts);
        opt_time = toc(opt_start);
        fprintf('\b\b\b\b: done in %g seconds\r', opt_time);
        
        % a posteriori checking
        Wf_opt = value(W_f);
        Wm_opt = value(W_m);
        Rus_opt = value(R_us);
        Rds_opt = value(R_ds);
        
        decided_vars = {Wf_opt, Wm_opt, Rus_opt, Rds_opt};
        
        % check constraint satisfaction
        failed = p.evaluate(ac, wind, t, decided_vars);
        failed_perc = failed / length(p.checked_constraints);
        
        % calculate R and rank W_s for every scenario 
        R = zeros(ac.N_G, N);
        ranks = zeros(N+1, 1);
        difference = zeros(N, 1);
        for i = 1:N
            for j = 1:ac.N_G
                k = ac.Gens(j);
                R(j, i) = trace(ac.Y_k(k) * (Wm_opt * wind.P_m(t,i))) - ac.C_w(k)*wind.P_m(t, i);
            end
            ranks(i) = svd_rank(Wf_opt + Wm_opt * wind.P_m(t,i));
            difference(i) = sum(R(:,i))+wind.P_m(t, i);
        end
        ranks(end) = svd_rank(Wf_opt);
        
        % calculate percentage of system matrices with rank 1
        rank_perc = sum(ranks==1)/(N+1)*100;
        
        % store results
        r = r + 1;
        results(r).formulation = 'P3';
        results(r).network = networks{n};
        results(r).scenarios = N;
        results(r).lambda = lambda;
        results(r).prep_time = prep_time;
        results(r).opt_time = opt_time;
        results(r).failed = failed_perc;
        results(r).mean_abs_diff = meanabs(difference);
        results(r).rank_perc = rank_perc;
        results(r).exit_code = diagnostics.problem;        
        results(r).decided_vars = decided_vars;
    end
end

%% TEST P3
p = formulate('P3*');
lambda = 0.5; % use fixed lambda for this

% Loop over different networks
for n = 1:length(networks)
    
    % load model
    ac = AC_model(networks{n});
    ac.set_WPG_bus(wind_bus(n));
    
    
    % define wind model with 24 timesteps and 20 percent WPG penetration
    wind = wind_model(ac, 24, 0.2);
    
    % define decision variables
    W_f = sdpvar(2*ac.N_b); 
    W_m = sdpvar(2*ac.N_b);
    R_us = sdpvar(ac.N_G, 1);
    R_ds = sdpvar(ac.N_G, 1);
    
    % loop over different scenario sizes
    for N = Ns
        
        fprintf('==================================================\nPreparing P3* %s with %i scenarios ...', networks{n}, N);

        % generate dummy scenarios
        wind.dummy(N);

        decision_vars = {W_f, W_m, R_us, R_ds};

        % prepare problem
        prep_start = tic;
        [Obj, Cons] = p.prepare(ac, wind, t, decision_vars, lambda);
        prep_time = toc(prep_start);
        fprintf('\b\b\b\b: done in %g seconds\n', prep_time);
        
        % optimize
        opts = sdpsettings('verbose', 0, 'solver', 'mosek');
        fprintf('Optimizing ...');
        opt_start = tic;
        diagnostics = optimize(Cons, Obj, opts);
        opt_time = toc(opt_start);
        fprintf('\b\b\b\b: done in %g seconds\r', opt_time);
        
        % a posteriori checking
        Wf_opt = value(W_f);
        Wm_opt = value(W_m);
        Rus_opt = value(R_us);
        Rds_opt = value(R_ds);
        
        decided_vars = {Wf_opt, Wm_opt, Rus_opt, Rds_opt};
        
        % check constraint satisfaction
        p_eval = formulate('P3'); % check against P3 constraints (psd for every scenario)
        failed = p_eval.evaluate(ac, wind, t, decided_vars);
        failed_perc = failed / length(p_eval.checked_constraints);
        
        % calculate R and rank W_s for every scenario 
        R = zeros(ac.N_G, N);
        ranks = zeros(N+1, 1);
        difference = zeros(N, 1);
        for i = 1:N
            for j = 1:ac.N_G
                k = ac.Gens(j);
                R(j, i) = trace(ac.Y_k(k) * (Wm_opt * wind.P_m(t,i))) - ac.C_w(k)*wind.P_m(t, i);
            end
            ranks(i) = svd_rank(Wf_opt + Wm_opt * wind.P_m(t,i));
            difference(i) = sum(R(:,i))+wind.P_m(t, i);
        end
        ranks(end) = svd_rank(Wf_opt);
        
        % calculate percentage of system matrices with rank 1
        rank_perc = sum(ranks==1)/(N+1)*100;
        
        % store results
        r = r + 1;
        results(r).formulation = 'P3*';
        results(r).network = networks{n};
        results(r).scenarios = N;
        results(r).lambda = lambda;
        results(r).prep_time = prep_time;
        results(r).opt_time = opt_time;
        results(r).failed = failed_perc;
        results(r).mean_abs_diff = meanabs(difference);
        results(r).rank_perc = rank_perc;
        results(r).exit_code = diagnostics.problem;        
        results(r).decided_vars = decided_vars;
    end
end