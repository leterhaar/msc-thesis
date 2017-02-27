%% runs the ADMM algorithm
clc


% load data
if not(exist('svd_rank', 'file'))
    addpath('../misc', '../networks', '../wind');
end
load('data/sdp_lib.mat')
% 
% slack = 1e-4;
% u_ = u_ + slack;
% l_ = l_ - slack;


% get dimensions from data
d = size(M_0, 1);
p = length(M_);
q = length(C_);
% M_0 = zeros(d);

it = struct();

% initiate variables
% Block 1
X_init = -ones(d);
it.X = X_init;
it.z_ = zeros(1,p);
it.z_0 = 0;

% Block 2
it.X_C_ = cell(1, q);
[it.X_C_{:}] = deal(X_init);
it.y_0 = 0;
it.y_ = zeros(1, p);
it.X_N_0 = X_init;
it.X_N_ = cell(p,1);
[it.X_N_{:}] = deal(X_init);

% Multipliers
it.lambda_z_ = zeros(1,p);
it.lambda_z_0 = 0;
it.Lambda_N_ = cell(1,p);
[it.Lambda_N_{:}] = deal(zeros(d));
it.Lambda_N_0 = zeros(d);
it.Lambda_C_ = cell(q,1);
[it.Lambda_C_{:}] = deal(zeros(d));

it.pres = nan;
it.energy = nan;
it.obj = nan;
it.obj2 = nan;
k = 1;

%% start iterations
max_its = 500;
prog = progress('Running ADMM', max_its-1);
while k < max_its
    mu = 1000;
    it(k+1).energy = 0;
    it(k+1).pres = 0; % primal residue
    
    %% Upate block 1
    % update X
    num = 0;
    den = 0;
    for r = 1:q
        num = num + C_{r} .* (it(k).X_C_{r} - (it(k).Lambda_C_{r}/mu));
        den = den + C_{r};
    end
    
    for s = 1:p
        num = num + N_{s} .* (it(k).X_N_{s} - (it(k).Lambda_N_{s}/mu));
        den = den + N_{s};
    end
    
    it(k+1).X = sparse_division(num, C, den);
    
    % update z_0
    it(k+1).z_0 = real(trace(M_0' * it(k).X_N_0) - ((it(k).lambda_z_0 + 1)/mu));

    % update z_s
    for s = 1:p
        it(k+1).z_(s) = max(min(real(trace(M_{s}' * it(k).X_N_{s}) ...
                                - (it(k).lambda_z_(s)/mu)), u_(s)), l_(s));
    end
    
    
    %% Update block 2
    
    % update X_C_r
    for r = 1:q
        X_indef = (it(k+1).X .* C_{r}) + (it(k).Lambda_C_{r}/mu);
        it(k+1).X_C_{r} = project_PSD(X_indef, C_{r});
%         it(k+1).energy = it(k+1).energy + mu*norm(it(k+1).X_C_{r} - it(k).X_C_{r}, 'fro')^2;
%         it(k+1).pres = it(k+1).pres + norm(it(k+1).X .* C_{r} - it(k+1).X_C_{r}, 'fro')^2;
    end
    
    % update y_0 and X_N_0
    it(k+1).y_0 = (it(k+1).z_0 + (it(k).lambda_z_0/mu) ...
                    - trace(M_0' * ((N_0 .* it(k+1).X) + ...
                    (it(k).Lambda_N_0/mu))))/(1+norm(M_0, 'fro')^2);
    it(k+1).X_N_0 = (N_0 .* it(k+1).X) + (it(k).Lambda_N_0/mu) + ...
                                                        (it(k+1).y_0 * M_0);
    it(k+1).energy = it(k+1).energy + mu*norm(it(k+1).X_N_0 - it(k).X_N_0, 'fro')^2;
    it(k+1).pres = it(k+1).pres + norm(it(k+1).X .* N_0 - it(k+1).X_N_0, 'fro')^2;

    
    % update y_s and X_N_s
    for s = 1:p
        it(k+1).y_(s) = (it(k+1).z_(s) + (it(k).lambda_z_(s)/mu) ...
                        - trace(M_{s}' * ((N_{s} .* it(k+1).X) + ...
                        (it(k).Lambda_N_{s}/mu))))/(1+norm(M_{s}, 'fro')^2);
        it(k+1).X_N_{s} = (N_{s} .* it(k+1).X) + (it(k).Lambda_N_{s}/mu) + ...
                                                     (it(k+1).y_(s) * M_{s});
%         it(k+1).energy = it(k+1).energy + mu*norm(it(k+1).X_N_{s} - it(k).X_N_{s}, 'fro')^2;
%         it(k+1).pres = it(k+1).pres + norm(it(k+1).X .* N_{s} - it(k+1).X_N_{s}, 'fro')^2;

    end
    
    %% Update multipliers
    
    % update Lambda_C_r
    for r = 1:q
        it(k+1).Lambda_C_{r} = it(k).Lambda_C_{r} + ...
                               mu * ((it(k+1).X .* C_{r}) - it(k+1).X_C_{r});
%         it(k+1).energy = it(k+1).energy + 1/mu * norm(it(k+1).Lambda_C_{r} - it(k).Lambda_C_{r}, 'fro')^2;
    end
    
    % update Lambda_N_s
    for s = 1:p
        it(k+1).Lambda_N_{s} = it(k).Lambda_N_{s} + ...
                               mu * ((it(k+1).X .* N_{s}) - it(k+1).X_N_{s});
%         it(k+1).energy = it(k+1).energy + 1/mu * norm(it(k+1).Lambda_N_{s} - it(k).Lambda_N_{s}, 'fro')^2;
    end
    
    % update Lambda_N_0
    it(k+1).Lambda_N_0 = it(k).Lambda_N_0 + ...
                               mu * ((it(k+1).X .* N_0) - it(k+1).X_N_0);
    it(k+1).energy = it(k+1).energy + 1/mu * norm(it(k+1).Lambda_N_0 - it(k).Lambda_N_0, 'fro')^2;
    
    % update lambda_z_s
    for s = 1:p
        it(k+1).lambda_z_(s) = real(it(k).lambda_z_(s) + ...
                               mu * (it(k+1).z_(s) - ...
                               trace(M_{s}' * it(k+1).X_N_{s})));
    end
%     it(k+1).energy = it(k+1).energy + 1/mu * norm(it(k+1).lambda_z_ - it(k).lambda_z_)^2;
    
    % update lambda_z_0
    it(k+1).lambda_z_0 = real(it(k).lambda_z_0 + ...
                               mu * (it(k+1).z_0 - ...
                               trace(M_0' * it(k+1).X_N_0)));
%     it(k+1).energy = it(k+1).energy + 1/mu * norm(it(k+1).lambda_z_0 - it(k).lambda_z_0)^2;
    
    % store objective
    it(k+1).obj = real(trace(M_0' * it(k+1).X_N_0));
    it(k+1).obj2 = real(trace(M_0' * it(k+1).X));
    
    % update iteration number
    k = k + 1;
    prog.ping();
    

    
    % check stopping criterion
%     if it(k).pres < 1e-10
%         prog.finish();
%         break;
%     end

    if k > 1
    % remove all items except X from memory to prevent crashing
        it(k-1).z_ = [];
        it(k-1).z_0 = [];

        it(k-1).X_C_ = [];
        it(k-1).y_0 = [];
        it(k-1).y_ = [];
        it(k-1).X_N_0 = [];
        it(k-1).X_N_ = [];

        % Multipliers
        it(k-1).lambda_z_ = [];
        it(k-1).lambda_z_0 = [];
        it(k-1).Lambda_N_ = [];
        it(k-1).Lambda_N_0 = [];
        it(k-1).Lambda_C_ = [];
    end
    
end
K = k;

%%
initfig('Residuals', 7);
hold off
semilogy([it.pres]); 
hold on;
semilogy([it.energy]);
grid on; box on;
legend('Primal', 'Energy');
xlabel('Iterations');
title('Residuals');

initfig('Objectives', 8);
plot(real([it.obj]));
plot(real([it.obj2]), '--')
plot(ones(k,1)*real(trace(M_0 * Xopt)));

legend('ADMM', 'ADMM - obj2', 'Centralized');
hold on; grid on; box on;

%% test solution against original problem

X = sdpvar(d, d, 'hermitian', 'complex');
Obj = trace(M_0' * X);
C = [];
for s = 1:length(l_)
    if isa(trace(X * M_{s}), 'sdpvar')
        C = [C; l_(s) - slack <= trace(M_{s}' * X) <= u_(s) + slack];
    end
end
for r = 1:length(C_)
    C = [C; X .* C_{r} >= 0];
end

assign(X, it(K).X);
if any(check(C) < -1e-4)
    fprintf('Infeasible (sum g<0 = %g)\n', sum(check(C(check(C) < -1e-4))));
else
    fprintf('Feasible\n');
end

fprintf('Rank X: %i\n', svd_rank(it(k).X, 1e1));
