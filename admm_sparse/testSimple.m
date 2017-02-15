% script to test whether the tree decomposition of the 14 bus network is
% the same as the full problem
% clear
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
t = 1; % timestep used for this demonstration (todo: add for everything)
N = 5;

% load network and wind models
ac = AC_model('case14a');
ac.set_WPG_bus(9);
wind = wind_model(ac, 24, 0.2);
wind.generate(N);

% add forecast as first scenario
wind.P_m = [zeros(24,1) wind.P_m];
wind.P_w = [wind.P_wf wind.P_w];


% optimization settings
ops = sdpsettings('solver', 'mosek');
%% define optimization variables and formulate objective

% define X \in H^{n*(N+1)}
X = sdpvar(ac.N_b*(N+1), ac.N_b*(N+1), 'hermitian', 'complex');

% define indices that correspond to the network states
select = cell(N+1,2);
for i = 1:N+1
    select{i,1} = ((i-1)*ac.N_b)+1:i*ac.N_b;
    select{i,2} = select{i,1};
end
% define R_us, R_ds and d_us, d_ds
R_us = sdpvar(ac.N_G,1);
R_ds = sdpvar(ac.N_G,1);
d_us = sdpvar(ac.N_G,1);
d_ds = sdpvar(ac.N_G,1);

%% define cost function
Obj = 0;

% real power generation in scenario 1 (= forecast)
for j = 1:ac.N_G
    k = ac.Gens(j);
    Obj = Obj + ac.c_qu(j)*((trace(X(select{1,:})*ac.Y_P(k))+ac.P_D(t,k))^2) ...
              + ac.c_li(j)*(trace(X(select{1,:})*ac.Y_P(k))+ac.P_D(t,k));
end

% add up and downspinning 
Obj = Obj + ac.c_us' * R_us + ac.c_ds' * R_ds;

%% define constraints (except psd)

C = [R_us >= 0, R_ds >= 0, sum(d_ds) == 1, sum(d_us) == 1];

% loop over scenario constraints
for i = 1:N+1
    % refbus angle constraints
    refbus_index = ac.refbus + (i-1)*ac.N_b;
    C = [C; imag(X(refbus_index, refbus_index)) == 0];
    
    for k = 1:ac.N_b
        % real power injection limits
        C = [C; (ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i) <= ...
                trace(X(select{i, :}) * ac.Y_P(k)) <= ...
                ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i)):...
                sprintf('Pinj | s%2i | b%2i', i, k)];

        % reactive power injection limits
        C = [C; (ac.Q_min(k) - ac.Q_D(t, k) <= ...
                trace(X(select{i, :}) * ac.Y_Q(k)) <= ...
                ac.Q_max(k) - ac.Q_D(t, k)):...
                sprintf('Qinj | s%2i | b%2i', i, k)];
        
        % voltage magnitude limits
        C = [C; ((ac.V_min(k))^2 <= ...
                trace(X(select{i, :}) * ac.M_k(k, 'single')) <= ...
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
                    trace((X(select{i, :}) - X(select{1, :})) * ac.Y_P(k)) ...
                    - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us(j)):...
                    sprintf('Rdus | s%2i | b%2i', i, k)];

            % relate W_s and W_f through d_ds and d_us
            C = [C; (trace((X(select{i, :}) - X(select{1, :})) * ac.Y_P(k)) ...
                    - ac.C_w(k)*wind.P_m(t, i) == ...
                    d_us(j) * max(0, -wind.P_m(t, i)) ...
                    - d_ds(j) * max(0, wind.P_m(t, i))):...
                    sprintf('Rbal | s%2i | b%2i', i, k)];        
        end
    end
end
%% solve problem with psd constraint on every scenario matrix variable
C_psdness_whole = [];

for i = 1:N+1
    
    C_psdness_whole = [C_psdness_whole; X(select{i,:}) >= 0];
end


tic
status = optimize([C_psdness_whole, C], Obj, ops);
toc
verify(not(status.problem), status.info);
verify(not(any(check(C) < -1e-6)), 'Infeasible solution!')

Xstar_whole = value(X);
Xstar_whole(isnan(Xstar_whole)) = 0; % remove nans
Rus_whole = value(R_us);
Rds_whole = value(R_ds);
dus_whole = value(d_us);
dds_whole = value(d_ds);
%% solve problem with psd constraint on trees

% manually enter bags for now, use algorithm from Madani & Lavaei later
bags = {[1,2,5],[2,4,5],[2,3,4],[4,5,9],[4,7,9],[7,8],[5,6,9],[6,9,13],...
        [9,13,14],[6,12,13],[6,9,11],[9,10,11]};

    
C_psdness = [];

for i = 1:N+1
    for j = 1:length(bags);
        
        bag = bags{j};
        bagsize = length(bag);
        selector = zeros(bagsize, (N+1)*ac.N_b);
        for k = 1:bagsize
            index = bag(k) + (i-1)*ac.N_b;   % select right entry
            selector(k, index) = 1;
        end
        C_psdness = [C_psdness; (selector * X * selector' >= 0):...
                        sprintf('PSD bag %2i | s %2i', j, i)];
    end
end

tic
status = optimize([C_psdness, C], Obj, ops);
toc
verify(not(status.problem), status.info);
verify(not(any(check(C) < -1e-6)), 'Infeasible solution!')
%%
Xstar_tree = value(X);

for i = 1:N+1
    verify(not(any(isnan(Xstar_tree(select{i,:})))), 'nans in tree submatrix for s%2i', i);
end
        
Xstar_tree(isnan(Xstar_tree)) = 0;
Rus_tree = value(R_us);
Rds_tree = value(R_ds);
dus_tree = value(d_us);
dds_tree = value(d_ds);
%% test rank of solutions and constraint satisfaction

verify(svd_rank(Xstar_whole) == N+1, 'whole matrix not rank %i', N+1);
verify(svd_rank(Xstar_tree) == N+1, 'tree matrix not rank %i', N+1);
verify(is_psd(Xstar_whole), 'whole matrix is not psd');
verify(is_psd(Xstar_tree), 'tree matrix is not psd');

for i = 1:N+1
    W_whole = Xstar_whole(select{i,:});
    verify(svd_rank(W_whole) == 1, sprintf('submatrix whole not rank 1 for s%2i', i));
    verify(is_psd(W_whole) == 1, sprintf('submatrix whole not psd for s%2i', i));
    W_tree = Xstar_tree(select{i,:});
    verify(svd_rank(W_tree) == 1, sprintf('submatrix tree not rank 1 for s%2i', i));
    verify(is_psd(W_tree) == 1, 'submatrix tree not psd for s%2i: %g', i, min(eig(W_tree)));
end

verify(all_close(Rus_tree, Rus_whole, 1e-6), 'Rus not close: %g', norm(Rus_tree-Rus_whole, 'inf'));
verify(all_close(Rds_tree, Rds_whole, 1e-6), 'Rds not close: %g', norm(Rds_tree-Rds_whole, 'inf'));
verify(all_close(dus_tree, dus_whole, 1e-6), 'dus not close: %g', norm(dus_tree-dds_whole, 'inf'));
verify(all_close(dds_tree, dds_whole, 1e-6), 'dds not close: %g', norm(dds_tree-dds_whole, 'inf'));

% compare PGs from full W and W with nans removed
for k = ac.Gens'
    PG1 = trace(ac.Y_P(k) * Xstar_whole(select{1, :}));
    PG2 = trace(ac.Y_P(k) * Xstar_tree(select{1, :}));
    [flag, distance] = all_close(PG1, PG2, 1e-6);
    verify(flag, 'PG%i not close: %g', k, distance)
end
