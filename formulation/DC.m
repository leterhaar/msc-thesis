addpath('../wind');
addpath('../misc');

figure(1);
set(1, 'name', 'Network');
dock

figure(2);
set(2, 'name', 'Wind');
dock

%%
dc = DC_model('case14');
dc.set_WPG_bus(9);
dc.draw_network();
dc.c_us(4) = 5;

epsilon = 5e-2;                         % violation parameter
zeta = 5*dc.N_G;                        % Helly-dimsension
beta = 1e-5;                            % confidence parameter

% determine number of scenarios to generate based on Eq (2-4)
N_orig = ceil(2/epsilon*(zeta-1+log(1/beta)));
wind_orig = wind_model(dc, 24, 0.2);
wind_orig.generate(N_orig);
%%
t = 5;
wind = copy(wind_orig);
wind.use_extremes(t);
N = 2;
wind_orig.plot(t);

P_G = sdpvar(dc.N_G, 1, 'full');       % real power production 
R_us = sdpvar(dc.N_G, 1, 'full');      % upspinning reserve requirement
R_ds = sdpvar(dc.N_G, 1, 'full');      % downspinning reserve requirement
d_us = sdpvar(dc.N_G, 1, 'full');      % upspinning distribution vector
d_ds = sdpvar(dc.N_G, 1, 'full');      % downspinning distribution vector
%% Define objective function

Obj = 0;

% loop over generators and add generator cost at time t
for k = 1:dc.N_G
    Obj = Obj + dc.c_qu(k) * (P_G(k))^2 + ...
                            dc.c_li(k) * P_G(k);
%       Obj = Obj + dc.c_us(k) * P_G(k);
end

% add reserve requirements costs
Obj = Obj + (dc.c_us' * R_us + dc.c_ds' * R_ds);
%% Define deterministic constraints
Cons = [];

% define deterministic power injection vector  
P_injf = (dc.C_G * P_G - dc.P_D(t, :)' + dc.C_w * wind.P_wf(t));

% power balance constraints
Cons = [Cons, ...
               ones(1, dc.N_b) * P_injf == 0];

% generator limits
Cons = [Cons, ...
               dc.P_Gmin <= P_G <= dc.P_Gmax];

% line flow limits
Cons = [Cons, ...
               - dc.P_fmax <=  ...
               dc.B_f * [dc.B_bustildeinv * P_injf(1:end-1); 0] ...
               <= dc.P_fmax];

% Non-negativity constraints for reserve requirements       
Cons = [Cons, R_us >= 0, R_ds >= 0];

Cons = [Cons, ones(1,dc.N_G)*d_ds == 1, ones(1,dc.N_G)*d_us == 1];

%% Define scenario constraints

% loop over scenarios
for i = 1:N

    % define reserve power
    R = d_us * max(0, -wind.P_m(t, i)) - d_ds * max(0, wind.P_m(t, i));

    % define scenario power injection vector
    P_injs = dc.C_G * (P_G + R) + dc.C_w * wind.P_w(t, i) - dc.P_D(t, :)';
    
%     % generator limits
    Cons = [Cons, ...
                   dc.P_Gmin <= P_G + R <= dc.P_Gmax];

%     % line flow limits
    Cons = [Cons, ...
                   - dc.P_fmax <=  ...
                   dc.B_f * [dc.B_bustildeinv * P_injs(1:end-1); 0]...
                   <= dc.P_fmax];
% 
%     % reserve requirements constraints
    Cons = [Cons, -R_ds <= R <= R_us];

end
%% the same as the three cells above
decision_vars = {P_G, R_us, R_ds, d_us, d_ds};

% p = formulate('DC');
% [Obj, Cons] = p.prepare(dc, wind, t, decision_vars);

%% Optimize
opt = sdpsettings('verbose', 0);
diagnostics = optimize(Cons, Obj, opt);

%% Evaluate and show
PG_opt = value(P_G);

Rus_opt = value(R_us);
Rds_opt = value(R_ds);

dus_opt = normV(value(d_us));
dds_opt = normV(value(d_ds));

% if only up or downspinning is required, replace the NaNs with zeros
if any(isnan(dus_opt))
    dus_opt = zeros(dc.N_G, 1);
end
if any(isnan(dds_opt))
    dds_opt = zeros(dc.N_G, 1);
end


% calculate R for every scenario
R = zeros(dc.N_G, N);
for i = 1:N
    for j = 1:dc.N_G
        k = dc.Gens(j);
        R(j, i) = d_us(j)*max(0, -wind.P_m(t, i)) - d_ds(j)*max(0, wind.P_m(t, i));
    end
end

% output ranks and table with results
decided_vars = {PG_opt, Rus_opt, Rds_opt, dus_opt, dds_opt};
format_result(dc, wind, t, decided_vars, R);

if diagnostics.problem ~= 0
    cprintf('red', '%s\n\n', diagnostics.info);
else
    fprintf('%s\n\n',diagnostics.info);
end