%% DC MODEL
% Class that constructs and stores all the values needed 
% for the DC model 
%
% Usage: ac = AC_model([modelname])
classdef AC_dry_model < handle
    properties
        model_name;      % name of the model
        N_b;             % number of buses [-]
        N_G;             % number of generators [-]
        N_l;             % number of branches / lines [-]
        N_w;             % number of wind buses
        P_max;           % Real injected power max
        P_min;           % Real injected power min
        Q_max;           % Reactive injected power max
        Q_min;           % Reactive injected power min
        V_min;           % Minimum bus voltage [p.u.]
        V_max;           % Maximum bus voltage [p.u.]
        S_max;           % Line flow limits
        P_lmmax;         % Real line flow limits
        Y;               % Complex nodal admittance matrix (N_b x N_b)
        P_D;             % Real power demanded [MW]
        Q_D;             % Reactive power demanded [MVAr]
        Ystar;           % Complex conjugate of Y
        from_to;         % relates topology of lines to line index
        Gens;            % list of generator bus indices
        C_w;             % Wind power connection matrix (N_b x N_w)
        C_G;             % Generator connection matrix (N_b x N_G)
        c_qu;            % quadratic generator cost vector [$/MWh]
        c_li;            % linear generator cost vector [$/MWh]
        c_us;            % upspinning cost vector [$/MWh]
        c_ds;            % downspinning cost vector [$/MWh]
        refbus;          % id of the reference bus
    end
   
    methods
        %% constructor function, calls the load function
        function obj = AC_dry_model(model_name)
            
            if nargin < 1
                % load case14 by default
                obj.model_name = 'case14';
            else
                obj.model_name = model_name;
            end
            
            filename = ['../networks/' obj.model_name '-saved.mat'];
            if exist(filename, 'file')
                saved_model = load(filename, 'obj');
                
                % loop over properties and assign
                for propname = fieldnames(saved_model.obj)'
                    obj.(propname{1}) = saved_model.obj.(propname{1});
                end
                
            else
                obj.make_model();
            end
            
            
        end
      
       
        function make_model(obj)
            %% loads the model and creates all the system matrices
            % Based on [Zimmerman2015] chapter 3
            mpc = loadcase(obj.model_name);
                        
            % determine number of branches, buses, generators and wind
            obj.N_l = size(mpc.branch, 1);
            obj.N_b = size(mpc.bus, 1);
            obj.N_G = size(mpc.gen, 1);
            obj.N_w = 0;   
           
            % extract bus voltage limits
            obj.V_max = mpc.bus(:, 12);          % Max bus voltage [p.u.]
            obj.V_min = mpc.bus(:, 13);          % Min bus voltage [p.u.]
            
            % store topology of lines and generators
            obj.from_to = mpc.branch(:, 1:2);
            obj.Gens = mpc.gen(:, 1);
            obj.refbus = find(mpc.bus(:, 2) == 3);
            
            for x = []
            % preallocate space for admittance vectors and switching matrices
            C_t = zeros(obj.N_l, obj.N_b);
            C_f = zeros(obj.N_l, obj.N_b);
            Y_ff = zeros(obj.N_l, 1);
            Y_tf = zeros(obj.N_l, 1);
            Y_ft = zeros(obj.N_l, 1);
            Y_tt = zeros(obj.N_l, 1);
            
            % loop through branches
            for i = 1:obj.N_l

                % extract branch data
                branch = mpc.branch(i, :);
                from = branch(1);                % from bus number
                to = branch(2);                  % to bus number
                r_s = branch(3);                 % series resistance [p.u.]
                x_s = branch(4);                 % series reactance [p.u.]
                b_c = branch(5);                 % series susceptance [p.u.]
                tau = branch(9);                 % turns ratio [x]
                theta = branch(10);              % phase shift [degrees]
                
                % NB discuss this with Vahab
                if tau == 0
                    tau = 1;
                end

                % construct from and to switching matrices
                from_i = find(mpc.bus(:,1) == from);
                to_i = find(mpc.bus(:,1) == to);
                C_t(i, from_i) = 1;
                C_f(i, to_i) = 1;
                
                % construct admittance vectors
                j = sqrt(-1);                    % imaginary unit
                y_s = 1/(r_s + (j * x_s));         % series admittance
                Y_ff(i) = (y_s + j * b_c / 2 ) * ( 1 / ( tau^2) );
                Y_ft(i) = -y_s * (1/(tau * exp(-j * theta)));
                Y_tf(i) = -y_s * (1/(tau * exp(j * theta)));
                Y_tt(i) = y_s + j * b_c / 2;
                
            end
            
            % make shunt admittance vector
            G_sh = mpc.bus(:, 5);                   % shunt conductance
            B_sh = mpc.bus(:, 6);                   % shunt susceptance 
            Y_sh = G_sh + j * B_sh;                 % shunt admittance
            
            % create nodal admittance matrix
            Y_f = diag(Y_ff) * C_f + diag(Y_ft) * C_t;
            Y_t = diag(Y_tf) * C_f + diag(Y_tt) * C_t;
            obj.Y = C_f' * Y_f + C_t' * Y_t + diag(Y_sh);
            obj.Ystar = conj(obj.Y);
            
            end
            
            % create switching matrix for generators C_g
            obj.C_G = zeros(obj.N_b, obj.N_G);
            for i = 1:obj.N_G
                bus_i = mpc.gen(i, 1);
                j = mpc.bus(:,1) == bus_i;
                obj.C_G(j, i) = 1;
            end
            
            % extract generator limits
            P_Gmax = mpc.gen(:, 9);          % Real upper limit [MV]
            P_Gmin = mpc.gen(:, 10);         % Real lower limit [MV]
            Q_Gmax = mpc.gen(:, 4);          % React upper limit [MVAr]
            Q_Gmin = mpc.gen(:, 5);          % React lower limit [MVAr]
            
            % convert to maximum bus injection
            baseMVA = 1;
            % baseMVA = 1;
            obj.P_max = obj.C_G * P_Gmax / baseMVA;
            obj.P_min = obj.C_G * P_Gmin / baseMVA;
            obj.Q_max = obj.C_G * Q_Gmax / baseMVA;
            obj.Q_min = obj.C_G * Q_Gmin / baseMVA;
            
            % demanded power
            obj.P_D = mpc.bus(:, 3) / baseMVA;
            obj.Q_D = mpc.bus(:, 4) / baseMVA;
            
            % distribute along profile
            load('load_profile.mat');
            obj.P_D = LoadProf * obj.P_D';
            obj.Q_D = LoadProf * obj.Q_D';
            
            obj.Y = makeYbus(mpc);
            obj.Ystar = obj.Y';
            
            % add generator costs if of the right format (quadratic)
            if isfield(mpc, 'gencost')
                is_of_type_2 = sum(mpc.gencost(:,1)) == obj.N_G*2;
                is_of_order_3 = sum(mpc.gencost(:,4)) == obj.N_G*3;
                if is_of_type_2 && is_of_order_3
                    obj.c_qu = mpc.gencost(:,5); % linear costs [$/MWh]
                    obj.c_li = mpc.gencost(:,4); % quadratic costs [$/MWh]
                end
            end
            
            % check if filling costs has worked, otherwise fill with random
            % variables
            if isempty(obj.c_qu)
                warning(['Could not load generator costs from MATPOWER' ...
                                    ' model, filled with random instead']);
                obj.c_li = rand(obj.N_G, 1);  % linear costs [$/MWh]
                obj.c_qu = rand(obj.N_G, 1);  % quadratic costs [$/MWh]
            end
           
            % determine up and downspinning costs based on steepest cost
            % curve
            obj.c_us = (obj.c_qu.*(P_Gmax.^2)+obj.c_li.*P_Gmax) ...
                                                        ./ P_Gmax;
                                                    
            % downspinning costs are the less than upspinning costs
            obj.c_ds = 0.1*obj.c_us;
            
            % define line flow limits (can be zero)
            obj.S_max = mpc.branch(:, 6) ./ baseMVA;
            
            obj.P_lmmax = ones(obj.N_l,1);
            
            % save the model
            save(['../networks/' obj.model_name '-saved.mat'], 'obj');
            disp(['Saved ' obj.model_name ' for future use']);
        
        end
        
        function Y_Qk = Y_Q(obj, k)
            %% returns the bus admittance Y_(Q;k)
            % these parameter matrices are defined in Madani2015
            % and are needed to construct the real and reactive power
            % injection
            %
            % Parameters
            % ----------
            % k : the index of the bus
            %
            % Returns
            % Y_{Q;k}
            
            % construct standard basis vector
            e_k = obj.e(k);
            j = sqrt(-1);
            Y_Qk = (1/(2*j))*(obj.Ystar*e_k*e_k' - e_k*e_k'*obj.Y);
        end
        
        function Y_Pk = Y_P(obj, k)
            %% returns the bus admittance Y_(P;k)
            % these parameter matrices are defined in Madani2015
            % and are needed to construct the real and reactive power
            % injection
            %
            % Parameters
            % k : the index of the bus
            %
            % Returns
            % Y_{P;k}
            
            % construct standard basis vector
            e_k = obj.e(k);
            Y_Pk = 0.5 * (obj.Ystar*e_k*e_k' + e_k*e_k'*obj.Y);
        end
        
        function Yk = Y_k(obj, k)
            %% Returns Y_k from Lavei 2012 paper
            e_k = obj.e(k);
            
            Y_k = e_k*e_k'*obj.Y;
            Yk = 0.5 * [real(Y_k + Y_k.') imag(Y_k.' - Y_k)
                        imag(Y_k - Y_k.') real(Y_k + Y_k.')];
        end
        
        function Ybark = Ybar_k(obj, k)
            %% Returns Ybar_k from Lavei 2012 paper
            e_k = obj.e(k);
            
            Y_k = e_k*e_k'*obj.Y;
            Ybark = (-1/2) * [imag(Y_k + Y_k.') real(Y_k - Y_k.')
                            real(Y_k.' - Y_k) imag(Y_k + Y_k.')];
        end
        
        function Ylm = Y_lm(obj, k)
            %% returns Y_lm from the k-th branch (from Lavei 2012)
            
            l = obj.from_to(k, 1);
            m = obj.from_to(k, 2);
            e_l = obj.e(l);
            e_m = obj.e(m);
            shunt = 0; % not available through MATPOWER
            
            y_lm = -obj.Y(l,m);
            Y_lm = (shunt + y_lm) * e_l * e_l' - (y_lm) * e_l * e_m';
            Ylm = 0.5 * [real(Y_lm + Y_lm.') imag(Y_lm.' - Y_lm)
                         imag(Y_lm - Y_lm.') real(Y_lm + Y_lm.')];
        end
        
        function Ybarlm = Ybar_lm(obj, k)
            %% returns Y_lm from the k-th branch (from Lavei 2012)
            
            l = obj.from_to(k, 1);
            m = obj.from_to(k, 2);
            e_l = obj.e(l);
            e_m = obj.e(m);
            shunt = 0; % not available through MATPOWER
            
            y_lm = -obj.Y(l,m);
            Y_lm = (shunt + y_lm) * e_l * e_l' - (y_lm) * e_l * e_m';
            Ybarlm = -0.5 * [imag(Y_lm + Y_lm.') real(Y_lm - Y_lm.')
                             real(Y_lm.' - Y_lm) imag(Y_lm + Y_lm.')];
        end
        
        function Mk = M_k(obj, k, flag)
            %% returns M_k from Lavei 2012
            assert(k <= obj.N_b, 'k cannot be larger than number of buses');
            e_k = obj.e(k);
            Mk = [e_k*e_k' zeros(obj.N_b)
                  zeros(obj.N_b) e_k*e_k'];
            % for the Madani version, only return the first block
            if nargin > 2
                Mk = Mk(1:obj.N_b, 1:obj.N_b);
            end
        end

        function set_WPG_bus(obj, WPG_bus_index)
            %% set_WPG_bus(bus_index)
            % Sets the WPG bus to a new index, updates C_w
            %
            % Parameters
            % ----------
            % bus_index : the index of the bus infeed (scalar or vector)
            
            % determine number of buses
            obj.N_w = length(WPG_bus_index);
            
            % fill connection matrix for wind C_w
            obj.C_w = zeros(obj.N_b, obj.N_w);
            for i = 1:length(WPG_bus_index)
                obj.C_w(WPG_bus_index(i), i) = 1;
            end
            
        end
        
        function e_k = e(obj, k)
            %% returns standard basis vector
            e_k = zeros(obj.N_b, 1);
            e_k(k) = 1;
        end
        
        function Ek = E_k(obj)
            %% returns the matrix needed to fix the refbus angle
            Ek = zeros(2*obj.N_b);
            refbus_im = obj.N_b + obj.refbus;
            Ek(refbus_im, refbus_im) = 1;
        end
        
        % draws the network layout
        function draw_network(obj)
            
            clf
            dock
            
            % load data
            mpc = loadcase(obj.model_name);
            
            % construct adjacency matrix
            adj_matrix = zeros(obj.N_b);
            for i = 1:obj.N_l
                from = find(mpc.bus(:,1) == mpc.branch(i, 1));
                to = find(mpc.bus(:,1) == mpc.branch(i, 2));
                adj_matrix(from,to) = 1;
                adj_matrix(to,from) = 1;
            end
            
            % define labels
            labels = textscan(num2str(mpc.bus(:,1)'),'   %s');
            
            % find generator indices
            gens = zeros(obj.N_G,1);
            for i = 1:obj.N_G
                gens(i) = find(obj.C_G(:,i));
            end
            
            % find wind buses
            winds = zeros(obj.N_w,1);
            for i = 1:obj.N_w
                winds = find(obj.C_w(:,i));
            end
            
            % find overlapping
            genwinds = gens(gens == winds);
            
            % plot the graph and highlight special nodes
            G = graph(adj_matrix, labels{1});

            clf
            h = G.plot();
            highlight(h, gens, ...
                    'Marker', 's', 'MarkerSize', 7, 'NodeColor', 'r');
            highlight(h, winds, ...
                    'NodeColor', 'g', 'Marker', 'p', 'MarkerSize', 7);
            highlight(h, genwinds, ...
                    'NodeColor', 'r', 'Marker', 'p', 'MarkerSize', 7);
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            title(['Schematic overview of ' obj.model_name ' network']);

        end
        
        function res = Y_(obj, k)
            % shortcut for Y_k(k)
            res = obj.Y_k(k);
        end
    end
end


