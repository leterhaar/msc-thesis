% DC SIMULATOR
% Class to bundle all the code regarding the validation of the
% probabilistic properties of a solution according to the DC model
%
% Usage :  
% DCsim = DC_simulator(P_G, d_ds, d_us, [[N_t], tol])
% DCsim.simulate(P_w, P_wf);


classdef DC_simulator < handle    
    properties
        m;              % An instance of DC_model
        P_G;            % The generator dispatch (N_G x N_t)
        d_us;           % The upspinning distribution (N_G x N_t)
        d_ds;           % The downspinning distribution (N_G x N_t)
        N_t;            % The number of time sttol
        P_w;            % The wind scenarios (N x N_t)
        V_balance;      % Violation for balancing constraints (1 x N_t)
        V_gen;          % Violation for generator limits (N_G x N_t)
        V_pf;           % Violation for power flow limits (N_l x N_t)
        V_total;        % Violation total (1 x N_t)
        N;              % Number of simulations done
        tol;            % Tolerance
    end
    
    
    methods
        
        function obj = DC_simulator(m, P_G, d_ds, d_us, N_t, tol)
        % DC_SIMULATOR creates the simulator
        % 
        % Parameters
        % m : Instance of 'DC_model'
        % P_G : Generator dispatch (N_G x N_t)
        % d_ds : Downspinning reserve  distribution (N_G x N_t)
        % d_us : Upspinning reserve  distribution (N_G x N_t)
        % N_t : # of timesteps. Default 24
        % tol : tolerance. Default 1e-5
        %
        % Returns
        % DC_simulator : instance of 'DC_simulator'
        
            if nargin < 6
                obj.tol = 1e-5;
            else
                obj.tol = tol;
            end
        
            if nargin < 5
                obj.N_t = 24;
            else
                obj.N_t = N_t;
            end
            
            % store in class
            obj.m = m;
            obj.P_G = P_G;
            obj.d_ds = d_ds;
            obj.d_us = d_us;
            
            obj.reset();
            
        end
        
        function simulate(obj, P_w, P_wf)
        % SIMULATE simulates the model using a wind realization and stores
        % the constraint violations in the class
        %
        % Parameters
        % P_w : Wind realization (N_t)
        % P_wf : Wind forecast (N_t)
        
            % wind power mismatch
            P_m = P_w - P_wf;
            
            % loop over time
            for t = 1:obj.N_t
                
                % define reserve
                R = obj.d_us(:, t) * max(0, -P_m(t)) ...
                    - obj.d_ds(:, t) * max(0, P_m(t));
                
                % define injection vector
                P_inj = obj.m.C_G * (obj.P_G(:, t) + R) + ...
                        obj.m.C_w * P_w(t) - obj.m.P_D(:, t);
                    
                % test power balance
                balance_violated = abs(ones(1, obj.m.N_b)*P_inj) > obj.tol;
                
                % test generator upper limits
                genup_violated = obj.P_G(:, t) + R > ...
                                 obj.m.P_Gmax + obj.tol;
               
                % test generator lower limits
                gendown_violated = obj.P_G(:, t) + R < ...
                                       obj.m.P_Gmin - obj.tol;
                
                % test power flow upper limits
                powerflows = obj.m.B_f * ...
                                [obj.m.B_bustildeinv * P_inj(1:end-1); 0];
                pf_violated = abs(powerflows) > obj.m.P_fmax + obj.tol;
                               
                % store all violations
                obj.V_balance(t) = obj.V_balance(t) + balance_violated;
                obj.V_gen(t) = obj.V_gen(t) +  ...
                             max(max(genup_violated, gendown_violated));
                obj.V_pf(t) = obj.V_pf(t) + max(pf_violated);
                obj.V_total(t) = obj.V_total(t) + ...
                        max([max(genup_violated), max(gendown_violated), ...
                            max(pf_violated), balance_violated]);
                
            end
            % increase simulation number
            obj.N = obj.N + 1;

        end
        
        function reset(obj)
        % RESET resets the simulation state and stored violations
        
            % set number of simulations to 0
            obj.N = 0;
            
            % (re-)initialize result matrices
            obj.V_balance = zeros(1, obj.N_t);
            obj.V_gen = zeros(1, obj.N_t);
            obj.V_pf = zeros(1, obj.N_t);
            obj.V_total = zeros(1, obj.N_t);
            
        end
    end
    
end