% clear settings;
% settings.epL = 1000;
% settings.line_prob = [2]';
% results = OPF_Solver('case3lesieutre',settings);

% results = OPF_Solver('case6ww');

% clear settings;
% settings.epB = 10;
% results = OPF_Solver('case9',settings);

results = OPF_Solver('case14');

% results = OPF_Solver('case24_ieee_rts');

% clear settings;
% settings.epB = 0.1;
% results = OPF_Solver('case30',settings);

% clear settings;
% settings.epB = 10;
% results = OPF_Solver('case39',settings);

% results = OPF_Solver('case57');

% clear settings;
% settings.epB = 10;
% results = OPF_Solver('case118',settings);

% clear settings;
% settings.epB = 10;
% results = OPF_Solver('case118cap_mod',settings);

% clear settings;
% settings.epB = 0.1;
% settings.epL = 100;
% results = OPF_Solver('case300',settings);
% settings.line_prob = results.line_prob;
% results = OPF_Solver('case300',settings);

% clear settings;
% settings.epB = 0.1;
% settings.epL = 100;
% settings.line_prob = [38; 402];
% results = OPF_Solver('case300',settings);

% clear settings;
% settings.epB = 3500;
% settings.epL = 3000;
% results = OPF_Solver('case2383wp',settings);
% settings.line_prob = results.line_prob;
% results = OPF_Solver('case2383wp',settings);

% clear settings;
% settings.epB = 3500;
% settings.epL = 3000;
% settings.line_prob = [100; 101; 102; 103; 104; 130; 134; 819; 821];
% results = OPF_Solver('case2383wp',settings);

% clear settings;
% settings.epB = 1500;
% results = OPF_Solver('case2736sp',settings);

% clear settings;
% settings.epB = 1000;
% results = OPF_Solver('case2737sop',settings);

% clear settings;
% settings.epB = 1000;
% results = OPF_Solver('case2746wop',settings);

% clear settings;
% settings.epB = 1000;
% results = OPF_Solver('case2746wp',settings);

% clear settings;
% settings.alpha = 1;
% settings.line_prob = 'all';
% settings.epL = 10000;
% settings.tol_feas = 1.5 * 10^(-5);
% results = OPF_Solver('case3012wp',settings);

% clear settings;
% settings.alpha = -1.5;
% settings.line_prob = 'all';
% settings.epL = 10000;
% settings.tol_feas = 1.5 * 10^(-5);
% results = OPF_Solver('case3120sp',settings);
