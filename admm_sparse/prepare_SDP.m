d = 10;
yalmip('clear');
slack = 0;
%% build system data matrices

M_0 = rhermsparse(d);
M_0 = 0.5 * (M_0 + M_0');

u_ = nan(0);
l_ = nan(0);
M_ = cell(0);
C_ = cell(0);
C = zeros(d);

for i = 1:10
    l_(end+1) = -rand;
    u_(end+1) = rand;
    M_{end+1} = rhermsparse(d);
end
for i = 1:10
    l_(end+1) = rand - slack ;
    u_(end+1) = l_(end) + slack;
    M_{end+1} = rhermsparse(d);
end
C_ = {ones(d)};
C = ones(d);

N_ = cell(length(M_), 1);
for i = 1:length(l_)
    assert(all(size(M_{i}) == [d d]));
    N_{i} = double(M_{i} ~= 0);
end
N_0 = double(M_0 ~= 0);



%% attempt to solve the program using YALMIP
X = sdpvar(d, d, 'hermitian', 'complex');
Obj = trace(M_0' * X);
Cons = [];
for s = 1:length(l_)
    if isa(trace(X * M_{s}), 'sdpvar')
        Cons = [Cons; l_(s) <= trace(M_{s} * X) <= u_(s)];
    end
end
for r = 1:length(C_)
    Cons = [Cons; X .* C_{r} >= 0];
end

ops = sdpsettings('verbose', 0, 'solver', 'mosek');
status = optimize(Cons, Obj, ops);
assert(not(status.problem), 'Primair probleem werkt niet: %s', status.info);
value(Obj)
Xopt_orig = value(X);
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
    Cons = [Cons; l_(s) <= z_(s) <= u_(s)];
end

for r = 1:q
    Cons = [Cons; X .* C_{r} == X_C_{r}];
    Cons = [Cons; X_C_{r} >= 0];
end

status = optimize(Cons, Obj, ops);
assert(not(status.problem), 'Equivalent probleem werkt niet: %s', status.info);
value(Obj)
[close, diff] = all_close(value(X), Xopt_orig, 1e-4)
Xopt = value(X);

%% 
save('data/sdp_simple', 'M_0', 'M_', 'l_', 'u_', 'N_', 'N_0', 'C_', 'C', 'Xopt', 'Xopt_orig');

