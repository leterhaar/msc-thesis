addpath('../networks');
addpath('../wind');
addpath('../formulation');

N = 10;
% create AC problem
ac = AC_model('case14a');
ac.set_WPG_bus(9);
ac_wind = wind_model(ac, 24, 0.2);
ac_wind.generate(N);

% create DC problem
dc = DC_model('case14a');
dc.set_WPG_bus(9);
wind = wind_model(dc, 24, 0.2);
wind.generate(N);
t_wind = randi(24);
x_sdp = sdpvar(5*dc.N_G, 1, 'full');
delta_sdp = sdpvar(1,3, 'full');
m = 4;

% create SVM problem
N_svm = 10;
d = 5;
svm = create_SVM(d,N_svm);
test_sequence = {   'equivalence_svm', ...
                    'equivalence_solvers', ...
                    'ACC_SVM', ...
                    'ACC_DC'};