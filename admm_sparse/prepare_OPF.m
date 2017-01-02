%% Add paths

if not(exist('AC_model', 'file'))
    addpath('../wind', '../misc', '../networks');
end

%% load models

N = 100;
t = 1;
N_w = 1;
N_lin = 2^N_w + 1;

ac = AC_model('case14a');
wind = wind_model(ac);
wind.generate(N);
wind.use_extremes_and_posneg(t);

% add first scenario as zero
wind.P_w = [wind.P_wf wind.P_w];
wind.P_m = [zeros(24,1) wind.P_m];

d = (N_lin)*ac.N_b + 4*ac.N_G; % dimensions of x
%% build system data matrices

% build M_0 for the objective
Y_sum = zeros(ac.N_b);
for j = 1:ac.N_G
    k = ac.Gens(j);
    Y_sum = Y_sum + ac.c_us(j) * ac.Y_P(k);
end
M_0 = blkdiag(Y_sum, zeros(ac.N_b*2), diag(ac.c_us), diag(ac.c_ds));

u_s = nan(0);
l_s = nan(0);
M = cell(0);

for i = 1:N_lin
    % refbus angle zero (don't know how to do this......)
    % ASK VAHAB
    
    for k = 1:ac.N_b
        % real power generation limits
        l_s(end+1) = ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i);
        M{end+1} = blkdiag(zeros(ac.N_b*(i-1)), ac.Y_P(k), zeros(ac.N_b*(N_lin-i) + 4*ac.N_G));
        u_s(end+1) =  ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t, i);
%         spy(M{end});
%         pause

        % reactive power injection limits
        l_s(end+1) = ac.Q_min(k) - ac.Q_D(t, k);
        M{end+1} = blkdiag(zeros(ac.N_b*(i-1)), ac.Y_Q(k), zeros(ac.N_b*(N_lin-i) + 4*ac.N_G));
        u_s(end+1) = ac.Q_max(k) - ac.Q_D(t, k);
%         spy(M{end})
%         pause

        l_s(end+1) = (ac.V_min(k))^2;
        M{end+1} = blkdiag(zeros(ac.N_b*(i-1)), ac.M_k(k, 1), zeros(ac.N_b*(N_lin-i) + 4*ac.N_G));
        u_s(end+1) = (ac.V_max(k))^2;
%         spy(M{end})
%         pause
    end
    
    if i > 1
        % reserve balancing constraints
        for j = 1:ac.N_G
                    % bus index
            k = ac.Gens(j);

            % Bound R between R_us and R_ds
            
            % -inf <= -[P^m_i]^- d^us_k - R^us_k <= 0
            l_s(end+1) = -inf;
            M{end+1} = blkdiag(zeros(N_lin*ac.N_b), ...
                       zeros(j-1), -wind.P_mneg(t, i), ...
                       zeros(2*ac.N_G-j), zeros(j-1), -1, zeros(2*ac.N_G-j));
            u_s(end+1) = 0;
%             spy(M{end})
%             pause
            
            % 0 <= -[P^m_i]^+ d^ds_k + R^ds_k <= +inf
            l_s(end+1) = 0;
            M{end+1} = blkdiag(zeros(N_lin*ac.N_b), ...
                       zeros(ac.N_G + j-1), -wind.P_mpos(t, i), ...
                       zeros(ac.N_G-j), zeros(ac.N_G + j-1), -1, zeros(ac.N_G-j));
            u_s(end+1) = inf;
%             spy(M{end})
%             pause 
    
            % Tr(Y_k W_i - W_0) + [P^m_i]^+ d^ds_k + [P^m_i]^- d^us_k ==
            % P^m_{i,k}
            l_s(end+1) = ac.C_w(k)*wind.P_mpos(t, i) + ...
                         ac.C_w(k)*wind.P_mneg(t, i);
            M{end+1} = blkdiag(-ac.Y_P(k), zeros(ac.N_b*(i-2)), ac.Y_P(k), ...
                               zeros(ac.N_b*(N_lin-i)), ...
                               zeros(j-1), wind.P_mneg(t, i), zeros(ac.N_G-j), ...
                               zeros(j-1), wind.P_mpos(t, i), zeros(3*ac.N_G-j));
            u_s(end+1) = ac.C_w(k)*wind.P_mpos(t, i) + ...
                          ac.C_w(k)*wind.P_mneg(t, i);
%             spy(M{end});
%             pause
        end
    end
    
end

% sum dus = 1
l_s(end+1) = 1;
M{end+1} = blkdiag(zeros(N_lin*ac.N_b), eye(ac.N_G), zeros(3*ac.N_G));
u_s(end+1) = 1;

% sum dds = 1
l_s(end+1) = 1;
M{end+1} = blkdiag(zeros(N_lin*ac.N_b), zeros(ac.N_G), eye(ac.N_G), zeros(2*ac.N_G));
u_s(end+1) = 1;

