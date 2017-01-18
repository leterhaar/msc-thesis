%% AC OPF
% OtH 29-9-16
addpath('../misc');
addpath('../wind');
addpath('../formulation');

%% Load models
ac = AC_model('case14a');
ac.set_WPG_bus(9);

N_t = 24;
wind = wind_model(ac, N_t, 0.2);

% define sample complexity
N = 4; 
% wind.generate(N);
% wind.use_forecast();
wind.dummy(N);
t = 8; % for now, do a loop over 24 hours later

%% Make initial guess by solving deterministic problem
W_f_det = sdpvar(2*ac.N_b); 

Obj_det = 0;
% P_G for every generator
for j = 1:ac.N_G
    k = ac.Gens(j);
    Obj_det = Obj_det + ac.c_us(j)*trace(ac.Y_k(k)*W_f_det);
end     
C_det = [];
k_ref = ac.N_b + ac.refbus;

% Deterministic constraints
for k = 1:ac.N_b
    % P_inj (1)
    C_det = [C_det, ac.P_min(k) ...
        <= trace(ac.Y_k(k)*W_f_det) + ac.P_D(t, k) - ac.C_w(k)*wind.P_wf(t) <= ...
           ac.P_max(k)];

    % Q_inj (2)
    C_det = [C_det, ac.Q_min(k) ...
        <= trace(ac.Ybar_k(k)*W_f_det) + ac.Q_D(t, k) <=...
           ac.Q_max(k)];

    % V_bus (3)
    C_det = [C_det, ac.V_min(k)^2 ...
        <= trace(ac.M_k(k)*W_f_det) <= ...
           ac.V_max(k)^2];
end

% PSD constraint on W_f 
C_det = [C_det, W_f_det >= 0];

% refbus constraint
C_det = [C_det, W_f_det(k_ref, k_ref) == 0];

% Optimize
opt = sdpsettings('verbose', 0, 'debug', 1, 'solver', 'mosek');
diagnostics = optimize(C_det, Obj_det, opt);
Wf_init = value(W_f_det);

%% Define real problem
W_f = sdpvar(2*ac.N_b); 
% W_f is a symmetric real valued matrix

W_s = cell(N, 1);
for i = 1:N
    W_s{i} = sdpvar(2*ac.N_b);
end

% put together in large matrix variable
X = blkdiag(W_f, W_s{:});

% preallocate the data matrices
data_matrices = zeros((N+1)*ac.N_b*2, (N+1)*2*ac.N_b, 3*(N+1)*ac.N_b+N*ac.N_G+N+2);

% maximum up and downspinning vectors
R_us = sdpvar(ac.N_G, 1);
R_ds = sdpvar(ac.N_G, 1);

% distribution vectors
d_us = sdpvar(ac.N_G, 1);
d_ds = sdpvar(ac.N_G, 1);

%% Define objective

% Y_obj
Yf = 0;
Ys = cell(N, 1);
[Ys{:}] = deal(zeros(2*ac.N_b));
lambda = 1/(2*N);
for j = 1:ac.N_G
    k = ac.Gens(j);
    % tilde Y f
    Yf = Yf + ac.c_us(j) * ac.Y_k(k) - lambda * N * ac.c_us(j) * ac.Y_k(k);
    for i = 1:N
        % tilde Y s 
        Ys{i} = Ys{i} + lambda * ac.c_us(j) * ac.Y_k(k);
    end
end
Y_obj = blkdiag(Yf, Ys{:});

data_matrices(:,:,1) = Y_obj;
Obj = trace(Y_obj * X);
%% Define constraints
C = [];

% define slack
slack = sdpvar(6*(N+1)*ac.N_b, 1, 'full');
C = [C, slack(:) >= 0];

% Deterministic constraints
for k = 1:ac.N_b
    
    slack_offset = 6*(k-1); % offset for slack variables
    dm_offset = 3*(k-1)+1;    % offset for data matrices
    
    % build tildeY_k, tilde Ybar_k and tilde M_k for W^f
    tildeY_k = blkdiag(ac.Y_k(k), zeros(N * 2 * ac.N_b));
    data_matrices(:,:,dm_offset+1) = tildeY_k;
    tildeYbar_k = blkdiag(ac.Ybar_k(k), zeros(N * 2 * ac.N_b));
    data_matrices(:,:,dm_offset+2) = tildeYbar_k;
    tildeM_k = blkdiag(ac.M_k(k), zeros(N*2*ac.N_b));
    data_matrices(:,:,dm_offset+3) = tildeM_k;
    
    % P_inj (1)
    C = [C, trace(tildeY_k * X) == ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_wf(t) + slack(slack_offset+1)];
    C = [C, trace(tildeY_k * X) == ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_wf(t) - slack(slack_offset+2)];
    
    % Q_inj (2)
    C = [C, trace(tildeYbar_k * X) == ac.Q_min(k) - ac.Q_D(t, k) + slack(slack_offset+3)];
    C = [C, trace(tildeYbar_k * X) == ac.Q_max(k) - ac.Q_D(t, k) - slack(slack_offset+4)];
    
    % V_bus (3)
    C = [C, trace(tildeM_k * X) == ac.V_min(k)^2 + slack(slack_offset+5)];
    C = [C, trace(tildeM_k * X) == ac.V_max(k)^2 - slack(slack_offset+6)];
    
end

%%
for i = 1:N
        
    for k = 1:ac.N_b
        
        slack_offset = (i * 6 * ac.N_b) + 6 * (k-1);    % slack offset
        dm_offset = (i*3*ac.N_b) + 3*(k-1) + 1;         % data matrices offset
        
        % build tildeY_k, tilde Ybar_k and tilde M_k for W^f
        tildeY_k = blkdiag(zeros(i*2*ac.N_b), ac.Y_k(k), zeros((N-i) * 2 * ac.N_b));
        data_matrices(:,:,dm_offset + 1) = tildeY_k;
        tildeYbar_k = blkdiag(zeros(i*2*ac.N_b), ac.Ybar_k(k), zeros((N-i) * 2 * ac.N_b));
        data_matrices(:,:,dm_offset + 2) = tildeYbar_k;
        tildeM_k = blkdiag(zeros(i*2*ac.N_b), ac.M_k(k), zeros((N-i) * 2 * ac.N_b));
        data_matrices(:,:,dm_offset + 3) = tildeM_k;
    
        % P_inj (1')
        C = [C, trace(tildeY_k * X) == ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i) + slack(slack_offset+1)];
        C = [C, trace(tildeY_k * X) == ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i) - slack(slack_offset+2)];

        % Q_inj (2')
        C = [C, trace(tildeYbar_k * X) == ac.Q_min(k) - ac.Q_D(t, k) + slack(slack_offset+3)];
        C = [C, trace(tildeYbar_k * X) == ac.Q_max(k) - ac.Q_D(t, k) - slack(slack_offset+4)];

        % V_bus (3')
        C = [C, trace(tildeM_k * X) == ac.V_min(k)^2 + slack(slack_offset+5)];
        C = [C, trace(tildeM_k * X) == ac.V_max(k)^2 - slack(slack_offset+6)];
        
    end
end

for i = 1:N

    for j = 1:ac.N_G
        
        dm_offset = 3*(N+1)*ac.N_b+1+(i-1)*ac.N_G;
        % bus index
        k = ac.Gens(j);
        
        tildeYdist = blkdiag(-ac.Y_k(k), zeros((i-1)*2*ac.N_b), ac.Y_k(k), zeros((N-i)*2*ac.N_b));
        data_matrices(:,:,dm_offset+j) = tildeYdist;
        % relate W_s and W_f through d_ds and d_us
        C = [C, trace(tildeYdist * X)  == ac.C_w(k)*wind.P_m(t, i) + ...
            d_us(j) * max(0, -wind.P_m(t, i)) - d_ds(j) * max(0, wind.P_m(t, i))];
    end
    
end

% refbus constraint
for i = 0:N % also include forecast
    dm_offset = 3*(N+1)*ac.N_b+N*ac.N_G+2;
    tildeEref = blkdiag(zeros(i*2*ac.N_b), ac.E_k, zeros((N-i)*2*ac.N_b));
    data_matrices(:,:,dm_offset+i) = tildeEref;
    C = [C, trace(tildeEref * X) == 0];
end

% PSD constraint on X
C = [C, X >= 0];
%% Optimize using MOSEK
opt = sdpsettings('verbose', 1, 'debug', 1, 'solver', 'mosek');
diagnostics = optimize(C, Obj, opt);

%% Preprocessing ADMM algorithm

% vectorize all data matrices
c = vec(Y_obj);
A = zeros(((N+1)*2*ac.N_b)^2, size(data_matrices, 3)-1);
for i = 1:size(data_matrices, 3)-1
    A(:, i) = vec(data_matrices(:,:,i+1));
end
A = A';
x = vec(X);

% create initial guess by using solution of deterministic problem
Ws_init = cell(N,1);
for i = 1:N
    Ws_init{i} = Wf_init;
end
x_init = vec(blkdiag(Wf_init, Ws_init{:}));

% create sparsity pattern of all data matrices
sparsity = sum(abs(data_matrices),3) > 0;
plot(graph(sparsity))


% chordal extension??? graph is not yet chordal and also not connected...
% HOW TO DO THIS?

% find cliques
cliques = maximalCliques(sparsity - diag(diag(sparsity))); % no self loops allowed for clique finding alogrithm
p = size(cliques,2);
H_ = cell(p,1);
x_ = cell(p,1);

% construct H_k and x_k for every clique
for k = 1:p
    clique = cliques(:,k);
    entries = find(clique);
    N_clique = length(entries);
    E_k = zeros(N_clique, (N+1)*2*ac.N_b);
    for i = 1:N_clique
        E_k(i, entries(i)) = 1;
    end
    H_{k} = kron(E_k, E_k);
    x_{k} = H_{k} * x_init;
end

% w/o chordal extension, slack and R cannot be included

%% Evaluate
Wf_opt = value(W_f);
Ws_opt = zeros(2*ac.N_b, 2*ac.N_b, N);
for i = 1:N
    Ws_opt(:,:,i) = value(W_s{i});
end

dus_opt = normV(value(d_us));
dds_opt = normV(value(d_ds));

% if only up or downspinning is required, replace the NaNs with zeros
if any(isnan(dus_opt))
    dus_opt = zeros(ac.N_G, 1);
end
if any(isnan(dds_opt))
    dds_opt = zeros(ac.N_G, 1);
end

% replace really small values with zeros
dds_opt(abs(dds_opt) < 1e-4) = 0;
dus_opt(abs(dus_opt) < 1e-4) = 0;

% calculate R for every scenario
R = zeros(ac.N_G, N);
P_G = zeros(ac.N_G, 1);
for i = 1:N
    for j = 1:ac.N_G
        k = ac.Gens(j);
        R(j, i) = trace(ac.Y_k(k) * (Ws_opt(:,:,i)-Wf_opt)) - ac.C_w(k)*wind.P_m(t, i);
        
        % calculate P_G
        if i == 1
            P_G(j) = trace(ac.Y_k(k) * Wf_opt) + ac.P_D(t, k) - ac.C_w(k)*wind.P_wf(k);
        end
    end
end

Rus_opt = max(R, [], 2);
Rds_opt = min(R, [], 2);

decided_vars = {Wf_opt, Ws_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};

% output ranks and table with results
format_result(ac, wind, t, decided_vars, R);

% output solver state 
fprintf('%s in %g seconds\n\n',diagnostics.info, diagnostics.solvertime);
%% check constraints
fprintf('Checking constraints...');
checked = check(C);
fprintf(repmat('\b', 1, 23));

if any(checked < -1e06)
    fprintf('Some constraints not satisfied...');
end