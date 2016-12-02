% should give errors if f is not sdpvar
try
    prox(1,2,3,4);
    error_thrown = false;
catch
    error_thrown = true;
end
assert(error_thrown, 'Should throw error');

% test for vector valued problems
x = sdpvar(2,1);
f = [5 3]*x;
z = [10; 1];
alpha = 0.5;
f_prox = f + 1/(2*alpha)*norm(x-z)^2;
opt_settings = sdpsettings('verbose', 0);

% without constraints
optimize([], f_prox, opt_settings);
next_x = value(x);
next_x2 = prox(f, alpha, z, x);
assert(all_close(next_x, next_x2));

% with constraints
cons = [x(1) <= 3, x(2) <= 0];
optimize(cons, f_prox, opt_settings);
next_x = value(x);
next_x2 = prox(f, alpha, z, x, cons);
assert(all_close(next_x, next_x2));

% matrix variables
p = read_sdpa('../sdplib/hinf1.dat-s');

f = trace(-p.C*p.X);
z = randd(size(p.X, 1),1);
Z = z*z';
f_prox = f + 1/(2*alpha)*norm(p.X-Z, 'fro')^2;

% with constraints
cons = [p.X >= 0];
for i = 1:length(p.As)
    cons = [cons, trace(p.As{i}*p.X) == p.bs(i)];
end
optimize(cons, f_prox, opt_settings);
next_x = value(p.X);
next_x2 = prox(f, alpha, Z, p.X, cons);
assert(all_close(next_x, next_x2));
