%% Add paths

if not(exist('AC_model', 'file'))
    addpath('../wind', '../misc', '../networks');
end
yalmip('clear');
clear
clc

%% load models

N = 100;
t = 1;
N_w = 1;
slack = 0;

ac = AC_model('case14a');
ac.set_WPG_bus(9);
wind = wind_model(ac);
wind.generate(N);
wind.use_extremes_and_posneg(t);
N = 2^N_w + 1;

% add first scenario as zero
wind.P_w = [wind.P_wf wind.P_w];
wind.P_m = [zeros(24,1) wind.P_m];

d = (N)*ac.N_b + 4*ac.N_G; % dimensions of x
%% build system data matrices

% build M_0 for the objective
Y_sum = zeros(ac.N_b);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Y_sum = Y_sum + ac.c_us(j) * ac.Y_P(k);
end
M_0 = blkdiag(Y_sum, zeros(ac.N_b*(N-1) + 2*ac.N_G), diag(ac.c_us), diag(ac.c_ds));

u_ = nan(0);
l_ = nan(0);
M_ = cell(0);
C_ = cell(0);
C = zeros(d);

bags = {[1,2,5],[2,4,5],[2,3,4],[4,5,9],[4,7,9],[7,8],[5,6,9],[6,9,13],...
        [9,13,14],[6,12,13],[6,9,11],[9,10,11]};


for i = 1:N
    % refbus angle zero (don't know how to do this......)
    % is this a problem....????
    
    for k = 1:ac.N_b
        % real power generation limits
        l_(end+1) = ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i);
        M_{end+1} = blkdiag(zeros(ac.N_b*(i-1)), ac.Y_P(k), zeros(ac.N_b*(N-i) + 4*ac.N_G));
        u_(end+1) =  ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i);
%         spy(M_{end});
%         pause

        % reactive power injection limits
        l_(end+1) = ac.Q_min(k) - ac.Q_D(t, k);
        M_{end+1} = blkdiag(zeros(ac.N_b*(i-1)), ac.Y_Q(k), zeros(ac.N_b*(N-i) + 4*ac.N_G));
        u_(end+1) = ac.Q_max(k) - ac.Q_D(t, k);
%         spy(M_{end}) 
%         pause

        l_(end+1) = (ac.V_min(k))^2;
        M_{end+1} = blkdiag(zeros(ac.N_b*(i-1)), ac.M_k(k, 1), zeros(ac.N_b*(N-i) + 4*ac.N_G));
        u_(end+1) = (ac.V_max(k))^2;
%         spy(M_{end})
%         pause
    end
    
    if i > 1
        % reserve balancing constraints
        for j = 1:ac.N_G
                    % bus index
            k = ac.Gens(j);

            % Bound R between R_us and R_ds
            
            % -1e6 <= -[P^m_i]^- d^us_k - R^us_k <= 0
            l_(end+1) = -1e3;
            M_{end+1} = blkdiag(zeros(N*ac.N_b), ...
                       zeros(j-1), -wind.P_mneg(t, i), ...
                       zeros(2*ac.N_G-j), zeros(j-1), -1, zeros(2*ac.N_G-j));
            u_(end+1) = 0;
%             spy(M_{end})
%             pause
            
            % 0 <= -[P^m_i]^+ d^ds_k + R^ds_k <= +1e6
            l_(end+1) = 0;
            M_{end+1} = blkdiag(zeros(N*ac.N_b), ...
                       zeros(ac.N_G + j-1), -wind.P_mpos(t, i), ...
                       zeros(ac.N_G-j), zeros(ac.N_G + j-1), 1, zeros(ac.N_G-j));
            u_(end+1) = 1e3;
%             spy(M_{end})
%             pause 
    
            % Tr(Y_k W_i - W_0) + [P^m_i]^+ d^ds_k + [P^m_i]^- d^us_k ==
            % P^m_{i,k}
            l_(end+1) = ac.C_w(k)*wind.P_mpos(t, i) + ...
                         ac.C_w(k)*wind.P_mneg(t, i);
            M_{end+1} = blkdiag(-ac.Y_P(k), zeros(ac.N_b*(i-2)), ac.Y_P(k), ...
                               zeros(ac.N_b*(N-i)), ...
                               zeros(j-1), -wind.P_mneg(t, i), zeros(ac.N_G-j), ...
                               zeros(j-1), -wind.P_mpos(t, i), zeros(3*ac.N_G-j));
            u_(end+1) = ac.C_w(k)*wind.P_mpos(t, i) + ...
                          ac.C_w(k)*wind.P_mneg(t, i);
%             spy(M_{end});
%             pause
        end
    end
    
    % build C
    for j = 1:length(bags)
        tildeC = zeros(ac.N_b);
        bag = bags{j};
        for b1 = bag
            for b2 = bag
                tildeC(b1, b2) = 1;
            end
        end
        C_{end+1} = blkdiag(zeros(ac.N_b*(i-1)), tildeC, ...
                           zeros(ac.N_b*(N-i) + 4*ac.N_G));
%         spy(C{end})
%         pause
        C = C + C_{end};
    end
%         C_{end+1} = blkdiag(zeros(ac.N_b*(i-1)), ones(ac.N_b), zeros(ac.N_b*(N-i) + 4*ac.N_G));
%         C = C + C_{end};
end

% sum dus = 1
l_(end+1) = 1;
M_{end+1} = blkdiag(zeros(N*ac.N_b), eye(ac.N_G), zeros(3*ac.N_G));
u_(end+1) = 1;

% sum dds = 1
l_(end+1) = 1;
M_{end+1} = blkdiag(zeros(N*ac.N_b), zeros(ac.N_G), eye(ac.N_G), zeros(2*ac.N_G));
u_(end+1) = 1;

N_ = cell(length(M_), 1);
for i = 1:length(l_)
    assert(all(size(M_{i}) == [d d]));
    N_{i} = M_{i} ~= 0;
end
N_0 = M_0 ~= 0;

% % add positivity constraint on R
% C = C + blkdiag(zeros(ac.N_b*N), eye(4*ac.N_G));
% % 
for j = 1:4*ac.N_G
    C_{end+1} = blkdiag(zeros(ac.N_b*N), zeros(j-1), 1, zeros(4*ac.N_G-j));
    C = C + C_{end};
end
C = double(C ~= 0);

save('data/14bus_RS', 'M_0', 'M_', 'l_', 'u_', 'N_', 'N_0', 'C_', 'C', 'ac', 'wind', 'slack');

%% attempt to solve the program using YALMIP
X = sdpvar(d, d, 'hermitian', 'complex');
Obj = trace(M_0' * X);
C = [];
for s = 1:length(l_)
    if isa(trace(X * M_{s}), 'sdpvar')
        C = [C; l_(s) <= trace(M_{s} * X) <= u_(s)];
    end
end
for r = 1:length(C_)
    C = [C; X .* C_{r} >= 0];
end

ops = sdpsettings('verbose', 0, 'solver', 'mosek', 'debug', 1);
status = optimize(C, Obj, ops);
verify(not(status.problem), 'Primair probleem werkt niet: %s', status.info);
value(Obj)
Xopt = value(X);
%% test in ADMM formulation
p = length(l_);
q = length(C_);
z_0 = sdpvar(1, 1);
z_ = sdpvar(p, 1, 'full', 'complex');
X_C_ = cell(q,1);
for r = 1:q
    X_C_{r} = sdpvar(d, d, 'hermitian', 'complex');
end
X_N_ = cell(p, 1);
for s = 1:p
    X_N_{s} = sdpvar(d,d,'hermitian','complex');
end
X_N_0 = sdpvar(d, d, 'hermitian', 'complex');

Obj = z_0;
C = [X .* N_0 == X_N_0];
C = [C; trace(M_0 * X_N_0) == z_0];
for s = 1:p
    C = [C; X .* N_{s} == X_N_{s}];
    C = [C; trace(M_{s} * X_N_{s}) == z_(s)];
    C = [C; l_(s) - slack <= z_(s) <= u_(s) + slack];
end

for r = 1:q
    C = [C; X .* C_{r} == X_C_{r}];
    C = [C; X_C_{r} >= 0];
end

status = optimize(C, Obj, ops);
verify(not(status.problem), 'Equivalent probleem werkt niet: %s', status.info);
value(Obj)
[close, diff] = all_close(zero_for_nan(value(X)), zero_for_nan(Xopt), 1e-4)