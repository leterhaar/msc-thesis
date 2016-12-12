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
 
test_sequence = {'equivalence_solvers', 'ACC'}; % 'AC_agent'