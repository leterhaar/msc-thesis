addpath('../networks');
addpath('../wind');
addpath('../formulation');

N = 10;
ac = AC_model('case9');
ac.set_WPG_bus(9);
wind = wind_model(ac, 24, 0.2);
wind.generate(N);


test_sequence = {'AC_agent'};