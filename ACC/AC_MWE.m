%% CREATE A RANDOM SDP

% define matrix variable
n = 10;
X = sdpvar(n);

% create m random constraints
m = 1000;
all_A = rand(n, m); % data matrices
all_b = rand(m,2); % data vectors

Cons = [];
Red_Cons = [];
for i = 1:m
    a = all_A(:,i);
    A = random_symmetric_matrix(n);
    b1 = all_b(i,1);
    b2 = -all_b(i, 2);
    Cons = [Cons; trace(A*X) <= b1; trace(A*X) >= b2];
    Red_Cons = [Red_Cons; trace(A*X) <= 2*b1; trace(A*X) >= 2*b2];
end

Cons = [Cons; X >= 0];
c = rand(n,1);
C = random_symmetric_matrix(n);
Obj = trace(C*X);

diagnostics = optimize([Red_Cons, Cons], Obj, sdpsettings('verbose', 0));
Xstar = value(X);
display(diagnostics.info);

diagnostics = optimize(Cons, Obj, sdpsettings('verbose', 0));
Xstar2 = value(X);
display(diagnostics.info);

if all_close(Xstar, Xstar2, 1e-3) display('Solution is the same'); end

%% Example from https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-251j-introduction-to-mathematical-programming-fall-2009/readings/MIT6_251JF09_SDP.pdf

A1 = [1 0 1; 0 3 7; 1 7 5];
A2 = [0 2 8; 2 6 0; 8 0 4];
C  = [1 2 3; 2 9 0; 3 0 7];
b1 = 11;
b2 = 19;

X = sdpvar(3);

Obj = trace(C*X);

Cons = [trace(A1*X) <= b1; trace(A2*X) <= b2; X >= 0];

diagnostics = optimize(Cons, Obj, sdpsettings('verbose', 0));
display(diagnostics.info);
Xstar = value(X);

residuals = check(Cons); % = 0.0000   60.1954   -0.0000

% formulate problem with only first constraint
Cons2 = [trace(A1*X) <= b1; X >= 0];

diagnostics = optimize(Cons2, Obj, sdpsettings('verbose', 0));
display(diagnostics.info);
Xstar2 = value(X);

if all_close(Xstar2, Xstar, 1e-3)
    fprintf('All close\n');
end

%% test with SDPLIB problem
tol = 1e-3;
[C, As, bs, X] = read_sdpa('../sdplib/hinf1.dat-s');
Obj = trace(C*X);
Cons = [X >= 0];
m = length(As);
slack = 5*ones(m,1);
for i = 1:m
    Cons = [Cons, trace(As{i}*X) <= bs(i) + slack(i), ...
                  trace(As{i}*X) >= bs(i) - slack(i)];
end

diagnostics = optimize(Cons, Obj, sdpsettings('verbose', 0));
display(diagnostics.info);
Xstar = value(X);

% optimize again, with only active constraints
Cons2 = [X >= 0];
for i = 1:m
    residual = trace(As{i}*Xstar) - bs(i) - slack(i);
    if residual > -tol && residual < tol
        Cons2 = [Cons2, trace(As{i}*X) <= bs(i) + slack(i)];
    end
    residual = trace(As{i}*Xstar) - bs(i) + slack(i);
    if residual > -tol && residual < tol
        Cons2 = [Cons2, trace(As{i}*X) >= bs(i) - slack(i)];
    end
end
residuals = check(Cons);
Cons3 = Cons(residuals > -tol & residuals < tol);
diagnostics = optimize(Cons3, Obj, sdpsettings('verbose', 0));
display(diagnostics.info);
Xstar2 = value(X);

assert(all_close(Xstar2, Xstar, tol), 'Not close');

