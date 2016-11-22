% create models
addpath('../wind');
addpath('../networks');
ac = AC_model('case14a');
ac.set_WPG_bus(9);
N = 50;
wind = wind_model(ac, 24, 0.2);
wind.dummy(N);
t = 1;

% create vars
x = {   sdpvar(2*ac.N_b), ...       Wf
        sdpvar(2*ac.N_b), ...       Wmus
        sdpvar(2*ac.N_b), ...       Wmds
        sdpvar(2*ac.N_G, 1)}; ...   Rus and Rds
random_vecs = rand(2*ac.N_b, 3);
random_x = cell(1,4);
for i = 1:3
    random_x{i} = random_vecs(:,i) * random_vecs(:,i)';
end
random_x{4} = rand(2*ac.N_G, 1);

% solve problem and store xstar
C = AC_cons_det(x, ac, wind, t);
for i = 1:N
    C = [C, AC_cons_scen(x, ac, wind.slice(i), t)];
end
Obj = AC_f(x, ac, wind, t);
optimize(C, Obj, sdpsettings('verbose', 0));
xstar = values_cell(x);

% define sequence of files to test
test_sequence = {...
    'AC_f', ...
    'AC_g', ...
    'AC_cons_det',...
    'AC_cons_scen',...
    'AC_check',...
    'AC_active'};
