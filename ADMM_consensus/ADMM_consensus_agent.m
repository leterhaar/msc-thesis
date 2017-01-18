classdef ADMM_consensus_agent < handle

    properties
        dv;          % decision vars
        it;          % iterations
        solver;      % optimizer object
        k;           % iteration number
        mu;          % stepsize
    end
    
    methods
        
        function a = ADMM_consensus_agent(ac, wind, t, mu)
        % ADMM_consensus_agent(ac, wind, t, mu)
        % initialize the agent, formulate optimization problem
        %
        % ARGUMENTS
        % =========
        % ac            : AC model
        % wind          : wind model with 1 scenario
        % t             : time for wind (TODO extend to multi hour)
        % mu            : fixed step size for ADMM
        % 
        % RETURNS
        % =======
        % agent         : a instance of ADMM_consensus_agent
            init_timer = tic;
            
            % create decision variables
            a.dv.Wf = sdpvar(2*ac.N_b);        % forecasted state
            a.dv.Ws = sdpvar(2*ac.N_b);        % scenario state
            a.dv.d = sdpvar(2*ac.N_G, 1);      % distribution vectors
            a.dv.R = sdpvar(2*ac.N_G, 1);      % reserve bounds
            
            % create solver parameters
            averages = sdpvar(4*ac.N_b^2 + 4*ac.N_G, 1);
            multipliers = sdpvar(4*ac.N_b^2 + 4*ac.N_G, 1);
            
            % create objective function
            obj = 0;
            for j = 1:ac.N_G
                k = ac.Gens(j);
                obj = obj ...
                        + ac.c_qu(j) * (trace(ac.Y_(k)*a.dv.Wf) + ...
                                      ac.P_D(t,k)-ac.C_w(k)*wind.P_wf(t))^2 ...
                        + ac.c_li(j) * (trace(ac.Y_(k)*a.dv.Wf) + ...
                                      ac.P_D(t,k)-ac.C_w(k)*wind.P_wf(t));
            end
            obj = obj + [ac.c_us' ac.c_ds'] * a.dv.R + ...
                  mu/2 * ([vec(a.dv.Wf); a.dv.d; a.dv.R] ...
                          - averages + (multipliers/mu))' * ...
                         ([vec(a.dv.Wf); a.dv.d; a.dv.R] ...
                          - averages + (multipliers/mu));
            
            % create constraints
            C = [a.dv.R >= 0; ...
                 sum(a.dv.d(1:ac.N_G)) == 1; ...
                 sum(a.dv.d(ac.N_G+1:2*ac.N_G)) == 1];
            
            % refbus angle constraints
            refbus_index = ac.refbus + ac.N_b;
            C = [C; a.dv.Wf(refbus_index, refbus_index) == 0];
            C = [C; a.dv.Ws(refbus_index, refbus_index) == 0];


            % psd constraints
            C = [C; a.dv.Wf >= 0; a.dv.Ws >= 0];

            for k = 1:ac.N_b
                % real power injection limits Wf
                C = [C; (ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_wf(t) <= ...
                        trace(a.dv.Wf * ac.Y_(k)) <= ...
                        ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_wf(t)):...
                        sprintf('Pinj Wf | b%2i', k)];

                % reactive power injection limits Wf
                C = [C; (ac.Q_min(k) - ac.Q_D(t, k) <= ...
                        trace(a.dv.Wf * ac.Ybar_k(k)) <= ...
                        ac.Q_max(k) - ac.Q_D(t, k)):...
                        sprintf('Qinj Wf | b%2i', k)];

                % voltage magnitude limits Wf
                C = [C; ((ac.V_min(k))^2 <= ...
                        trace(a.dv.Wf * ac.M_k(k)) <= ...
                        (ac.V_max(k))^2):...
                        sprintf('Vbus Wf | b%2i', k)];
                    
                % real power injection limits Ws
                C = [C; (ac.P_min(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t) <= ...
                        trace(a.dv.Ws * ac.Y_(k)) <= ...
                        ac.P_max(k) - ac.P_D(t, k) + ac.C_w(k)*wind.P_w(t)):...
                        sprintf('Pinj Ws | b%2i', k)];

                % reactive power injection limits Ws
                C = [C; (ac.Q_min(k) - ac.Q_D(t, k) <= ...
                        trace(a.dv.Ws * ac.Ybar_k(k)) <= ...
                        ac.Q_max(k) - ac.Q_D(t, k)):...
                        sprintf('Qinj Ws | b%2i', k)];

                % voltage magnitude limits Ws
                C = [C; ((ac.V_min(k))^2 <= ...
                        trace(a.dv.Ws * ac.M_k(k)) <= ...
                        (ac.V_max(k))^2):...
                        sprintf('Vbus Ws | b%2i', k)];
            end

            % reserve balancing constraints
            for j = 1:ac.N_G
                % bus index
                k = ac.Gens(j);

                % Bound R between R_us and a.dv.R
                C = [C; (-a.dv.R(j + ac.N_G) <= ...
                        trace((a.dv.Wf - a.dv.Ws) * ac.Y_(k)) ...
                        - ac.C_w(k)*wind.P_m(t) <= ...
                        a.dv.R(j)):...
                        sprintf('Rdus | b%2i', k)];

                % relate W_s and W_f through d_ds and d_us
                C = [C; (trace((a.dv.Wf - a.dv.Ws) * ac.Y_(k)) ...
                        - ac.C_w(k)*wind.P_m(t) == ...
                        - a.dv.d(j) * min(0, wind.P_m(t)) ...
                        - a.dv.d(j+ac.N_G) * max(0, wind.P_m(t))):...
                        sprintf('Rbal | b%2i', k)];        
            end
            
            % define solver
            a.solver = optimizer(C, obj, sdpsettings('solver', 'mosek'), ...
                                 {averages, multipliers}, {a.dv.Wf, ...
                                                           a.dv.Ws, ...
                                                           a.dv.d, ...
                                                           a.dv.R});
                                                            
            % define iterations and define initialization
            x0 = [ones(4*ac.N_b^2,1); zeros(4*ac.N_G,1);];
            a.it.lambda = zeros(4*ac.N_b^2 + 4*ac.N_G, 1);
            a.it.x = x0;
            a.it.time = toc(init_timer);
            a.it(2).x = x0;
            
            a.k = 2;
            a.mu = mu;
            
        end
        
        function x = update(a, averages)
        % take a step in the ADMM algorithm
        % 
        % ARGUMENTS
        % =========
        % averages      : vector with averages
        %
        % RETURNS
        % x_i           : vector with new local estimate of x
        
            step_timer = tic;
            
            % store residual
            a.it(a.k).residual = norm(a.it(a.k).x - averages)^2;

            % update multipliers
            a.it(a.k).lambda = a.it(a.k-1).lambda ...
                               + a.mu * (a.it(a.k).x - averages);

            % update x
            [X, problem, msg] = a.solver(averages, a.it(a.k).lambda);
            verify(not(problem), msg{:});

            % store in vector and broadcast
            x = [vec(X{1}); X{3}; X{4}];

            a.it(a.k+1).x = x;
            a.it(a.k+1).time = toc(step_timer);
            
            a.k = a.k + 1;
        end
        
    end
end