%% AC OPF (deterministic)
% OtH 31-8-16
addpath('../misc');

%% Load model
m = AC_model('case14');

% use simpler Y from Lavei2012
% m.override_Y();
N_t = 24;

%% Define optimization variable
W = sdpvar(m.N_b, m.N_b, 'hermitian', 'complex'); 
% W is an Hermitian matrix

%% Define objective function
Obj = FrobInProd(W, (m.Y + m.Ystar)/2);

%% Define constraints

Cons = [];

I = eye(m.N_b);

% loop over buses
for k = 1:m.N_b
    
    % bus voltage limits
    e_k = I(:, k);
    Cons = [Cons; m.V_min(k)^2 <= FrobInProd(W, e_k*e_k') <= m.V_max(k)^2];
    
    % reactive power injection limits
    Cons = [Cons; m.Q_min(k) <= FrobInProd(W, m.Y_Q(k)) + m.Q_D(k) <= ...
                                                m.Q_max(k) ];
                                          
    % real power injection limits
    Cons = [Cons; m.P_min(k) <= FrobInProd(W, m.Y_P(k)) + m.P_D(k) <= ...
                                                m.P_max(k)];
                                                                                
end

% PSD constraint on W
Cons = [Cons; W >= 0];

%% Optimize
opt_settings = sdpsettings();
optimize(Cons, Obj, opt_settings)

% runs into numerical problems...

%% Extract generator dispatch 
Wopt = value(W);
P_G = zeros(m.N_b,1);
Q_G = zeros(m.N_b,1);

for k = 1:m.N_b
    Q_G(k) = FrobInProd(Wopt, m.Y_Q(k)) + m.Q_D(k);
    P_G(k) = FrobInProd(Wopt, m.Y_P(k)) + m.P_D(k);
end
P_G(P_G < 1e-10) = 0;
Q_G(Q_G < 1e-10) = 0;

% discard the imaginary parts (very small)
assert(sum(imag(P_G)) < 1e-10 && sum(imag(Q_G)) < 1e-10);
P_G = real(P_G);
Q_G = real(Q_G);

Qg = m.C_G' * Q_G;
Pg = m.C_G' * P_G;

display(sprintf('Total reactive power generated: %f\tdemand: %f',...
    sum(Qg), sum(m.Q_D)));
display(sprintf('Total     real power generated: %f\tdemand: %f', ...
    sum(Pg), sum(m.P_D)));
