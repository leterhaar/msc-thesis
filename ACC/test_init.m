% create models
addpath('../wind');
addpath('../formulation');
ac = AC_model('case9');
ac.set_WPG_bus(9);
wind = wind_model(ac, 24, 0.2);
wind.generate(3);

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

test_sequence = {...
    'AC_f', ...
    'AC_g', ...
    'AC_cons_det',...
    'AC_cons_scen',...
... 'AC_agent', 
    'AC_active'};
