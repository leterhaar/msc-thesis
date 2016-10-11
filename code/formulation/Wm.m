%% AC OPF
% OtH 6-9-16
addpath('../misc');
addpath('../wind');
clear all
%% Load models
ac = AC_model('case14');
ac.set_WPG_bus(4);
ac.draw_network
%%
wind = wind_model(ac, 24, 0.2);

% define sample complexity
N = 2; 
wind.generate(N);

%% Define problem
Wf = sdpvar(2*ac.N_b); 
% Wf is a symmetric real valued matrix
Wm = sdpvar(2*ac.N_b);
% Wm is a symmetric real valued matrix      
Rup = sdpvar(ac.N_G, 1);
Rdown = sdpvar(ac.N_G, 1);

problem = formulation('P4');

t = 10;
wind.plot(t);


[Obj, Cons] = problem.prepare(ac, Wf, Wm, wind.P_wf, wind.P_w, Rup, Rdown, t);


%% Optimize
display(optimize(Cons, Obj));

%% Evaluate
Wf_opt = value(Wf);
Wm_opt = value(Wm);
Rup_opt = value(Rup);
Rdown_opt = value(Rdown);

R = zeros(ac.N_b, N);
for i = 1:N
    for k = 1:ac.N_b
        R(k, i) = trace(ac.Y_k(k) * Wm_opt) - ac.C_w(k)*wind.P_m(t, i);
    end
end
display([ac.Gens R(ac.Gens, :)])

fprintf('Rank Wf: %i\n', svd_rank(Wf_opt));
fprintf('Rank Wm: %i\n\n', svd_rank(Wm_opt));

problem.evaluate(ac, Wf_opt, Wm_opt, wind.P_wf, wind.P_w, Rup_opt, Rdown_opt, t);