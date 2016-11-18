% initialize agent
agent = AC_agent(ac, wind, 1, 1, 3);

Ncons = 6*ac.N_b + 2*ac.N_G + 1;
% check size of constraints
assert(length(agent.C_1_params) == 3*Ncons, 'not right size');