%% Add paths

if not(exist('AC_model', 'file'))
    addpath('../wind', '../misc', '../networks');
end

%% load models
t = 1;
slack = 1e-6;
ac = AC_model('case14');

d = ac.N_b; % dimensions of x
%% build system data matrices

% build M_0 for the objective
Y_sum = zeros(ac.N_b);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Y_sum = Y_sum + ac.c_us(j) * ac.Y_P(k);
end
M_0 = Y_sum;
% M_0 = (ac.Y + ac.Ystar) ./ 2;
assert(all(all(M_0 == M_0')))

u_ = nan(0);
l_ = nan(0);
M_ = cell(0);
C_ = cell(0);
C = zeros(d);

bags = { [1,2,5], [2,3,4], [2,4,5], [7,8], [4,7,9], [6,12,13], [4,5,9], [5,6,9], [9,10,11], [6,9,11], [6,9,13], [9,13,14], [13,14], [14]};


    
for k = 1:ac.N_b
    % real power generation limits
    l_(end+1) = ac.P_min(k) - ac.P_D(t, k);
    M_{end+1} = ac.Y_P(k);
    u_(end+1) =  ac.P_max(k) - ac.P_D(t, k);
%         spy(M{end});
%         pause

    % reactive power injection limits
    l_(end+1) = ac.Q_min(k) - ac.Q_D(t, k);
    M_{end+1} = ac.Y_Q(k);
    u_(end+1) = ac.Q_max(k) - ac.Q_D(t, k);
%         spy(M{end})
%         pause

    l_(end+1) = (ac.V_min(k))^2;
    M_{end+1} = ac.M_k(k, 1);
    u_(end+1) = (ac.V_max(k))^2;
%         spy(M{end})
%         pause
end

p = length(l_);

% build C
for j = 1:length(bags)
    tildeC = zeros(ac.N_b);
    bag = bags{j};
    for b1 = bag
        for b2 = bag
            tildeC(b1, b2) = 1;
        end
    end
    C_{end+1} = tildeC;
%         spy(C{end})
%         pause
    C = C + C_{end};
end
C = double(C ~= 0);
% C_ = {ones(d)};
% C = ones(d);

q = length(C_);

N_ = cell(length(M_), 1);
for i = 1:length(l_)
    assert(all(size(M_{i}) == [d d]));
    N_{i} = double(M_{i} ~= 0);
end
N_0 = double(M_0 ~= 0);




%% attempt to solve the program using YALMIP
X = sdpvar(d, d, 'hermitian', 'complex');
Obj = trace(M_0 * X);
Cons = [];
for s = 1:length(l_)
    if isa(trace(X * M_{s}), 'sdpvar')
        Cons = [Cons; l_(s) - slack <= trace(M_{s} * X) <= u_(s) + slack];
    end
end
for r = 1:length(C_)
    Cons = [Cons; X .* C_{r} >= 0];
end

ops = sdpsettings('verbose', 0, 'solver', 'mosek');
status = optimize(Cons, Obj, ops);
assert(not(status.problem), 'Primair probleem werkt niet: %s', status.info);
value(Obj)
Xopt_orig = (value(X));

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
Cons = [X .* N_0 == X_N_0];
Cons = [Cons; trace(M_0 * X_N_0) == z_0];
for s = 1:p
    Cons = [Cons; X .* N_{s} == X_N_{s}];
    Cons = [Cons; trace(M_{s} * X_N_{s}) == z_(s)];
    Cons = [Cons; l_(s) - slack <= z_(s) <= u_(s) + slack];
end

for r = 1:q
    Cons = [Cons; X .* C_{r} == X_C_{r}];
    Cons = [Cons; X_C_{r} >= 0];
end

status = optimize(Cons, Obj, ops);
assert(not(status.problem), 'Equivalent probleem werkt niet: %s', status.info);
value(Obj)
Xopt = zero_for_nan(value(X));

[close, diff] = all_close(Xopt, Xopt_orig, 1e-4)

save('data/14bus', 'M_0', 'M_', 'l_', 'u_', 'N_', 'N_0', 'C_', 'C', 'ac', 'Xopt', 'Xopt_orig');


%% Checking again
% Quickly check if the normal OPF problem is indeed feasible and 
% has the same solution

