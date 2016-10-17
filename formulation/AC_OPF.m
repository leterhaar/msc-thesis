%% AC OPF (deterministic)
% OtH 31-8-16
addpath('../misc');

%% Load model
ac = AC_model('case14');

%% Define optimization variable
W = sdpvar(2*ac.N_b); 
% W is an 2N_b * 2N_b symmetric real valued matrix


%% Define objective function
Obj = 0;
for k = ac.Gens'
    
    % have not used a cost vector, just real generated power
    Obj = Obj + trace(ac.Y_k(k)*W) + ac.P_D(k);
end

%% Define constraints

Cons = [];

% reference bus imaginary part should be fixed to zero, shift with N_b
% since imaginary part comes second
Cons = [Cons, W(ac.N_b+ac.refbus, ac.N_b+ac.refbus) == 0];

% PSD constraint on W
Cons = [Cons, W >= 0];

% loop over buses
for k = 1:ac.N_b
    
    % bus voltage limits
    Cons = [Cons, (ac.V_min(k)^2 <= trace(ac.M_k(k)*W) <= ac.V_max(k)^2)];
    
    % real power injection limits
    Cons = [Cons, (ac.P_min(k) - ac.P_D(k) <= trace(ac.Y_k(k)*W) <= ...
                                                ac.P_max(k) - ac.P_D(k))];
    % reactive power injection limits
    Cons = [Cons, (ac.Q_min(k) - ac.Q_D(k) <= trace(ac.Ybar_k(k)*W) <= ...
                                                ac.Q_max(k) - ac.Q_D(k))];                                                                                   
end

%% Optimize
display(optimize(Cons, Obj));

% Extract generator dispatch 
Wopt = value(W);


%% Check all constraints 
validator = AC_validator(ac);

validator.OPF(Wopt);
