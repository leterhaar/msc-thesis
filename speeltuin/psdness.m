n = 50;
N = 100;
a = randn(n,1);
addpath('../misc');

% create rank 1 psd matrix
A = a*a';

% define second matrix
B = sdpvar(n);

% define scaling
s = linspace(-10, 10, N);
s_min = min(s);
s_max = max(s);

% define PSD constraints and optimize
Cons = [];
for i = 1:N
    Cons = [Cons, A+s(i)*B >= 0];
end
optimize(Cons, trace(B));
B_all = value(B);

% only define extremes and optimize
Cons = [A+s_min*B >= 0, A+s_max*B >= 0];
optimize(Cons, trace(B));
B_extreme = value(B);

% check the result is the same
assert(all_close(B_extreme, B_all));