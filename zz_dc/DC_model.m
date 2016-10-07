%% DC MODEL
% Class that constructs and stores all the values needed 
% for the DC model
%
% Usage: dcmodel = DC_model([modelname])

classdef DC_model < handle
    properties
        model_name;      % name of the model
        N_b;             % number of buses [-]
        N_G;             % number of generators [-]
        N_l;             % number of branches / lines [-]
        N_w;             % number of wind buses
        Gens;            % set of generator bus indices
        B_bus;           % Bus susceptance
        B_bustildeinv;   % Inverted bus susceptance without last row
        B_f;             % Branch subsceptance
        P_busshift;      % Real power due to angle shift at bus [MW]
        P_fshift;        % Real power due to angle shift at branch [MW]
        P_Gmin;          % Minimal generator power [MW]
        P_Gmax;          % Maximal generator power [MW]
        G_sh;            % Shunt conductance [MW] (demanded at V = 1 p.u.)
        P_D;             % Demanded real power [MW]
        P_fmax;          % Maximum line rating [MW] (!!!)
        C_G;             % Generator connection matrix (N_b x N_G)
        C_w;             % Wind power connection matrix (N_b x N_w)
        c_qu;            % quadratic generator cost vector [$/MWh]
        c_li;            % linear generator cost vector [$/MWh]
        c_us;            % upspinning cost vector [$/MWh]
        c_ds;            % downspinning cost vector [$/MWh]
    end
   
    methods
        % constructor function, calls the load function
        function obj = DC_model(modelname)
            if nargin > 0
               obj.load_model(modelname);
               
            else
                % load case14 by default
                obj.load_model('case14');
            end
        end
      
        % loads the model and creates all the system matrices    
        function load_model(obj, modelname)
            obj.model_name = modelname;
            mpc = loadcase(modelname);
            
            % determine number of branches and buses
            obj.N_l = size(mpc.branch, 1);
            obj.N_b = size(mpc.bus, 1);
            obj.N_G = size(mpc.gen, 1);
            obj.N_w = 0;   % hardcoded for now, can be changed
            
            % find reference bus index
            refbus = find(mpc.bus(:, 2) == 3); 
            obj.Gens = mpc.gen(:, 1);

            
            % reorder buses to let refbus be the last row
            mpc.bus = mpc.bus([find(1:obj.N_b ~= refbus) ...
                                refbus], :);
            
            % extract bus, branch and gen data
            obj.G_sh = mpc.bus(:, 5);            % Shunt conductance
            obj.P_Gmax = mpc.gen(:, 9);          % Generator max limit [MV]
            obj.P_Gmin = mpc.gen(:, 10);         % Generator min limit [MV]

            % preallocate space for system matrices and vectors
            B_ff = zeros(obj.N_l, 1);            % branch susceptance
            C_t = zeros(obj.N_l, obj.N_b);       % connection matrix to
            C_f = zeros(obj.N_l, obj.N_b);       % connection matrix from
            obj.P_fshift = zeros(obj.N_l, 1);    % power due to angle shift

            % loop through branches
            for i = 1:obj.N_l

                % extract branch data
                branch = mpc.branch(i, :);
                tau = branch(9);                 % turns ratio [x]
                x = branch(4);                   % series reactance [p.u.]
                theta = branch(10);              % phase shift [degrees]
                from = branch(1);                % from bus number
                to = branch(2);                  % to bus number

                % if turn ratio is zero... 
                if tau == 0
                    % ... the susceptance is approximated by 
                    % the inverse of the reactance
                    B_ff(i) = 1/x;
                else
                    % else, the susceptance is approximated by the inverse 
                    % of the reactance multiplied by the turn ratio
                    B_ff(i) = 1/(x*tau);
                end

                obj.P_fshift(i) = -theta * B_ff(i);

                % construct from and to switching matrices
                from_i = find(mpc.bus(:,1) == from);
                to_i = find(mpc.bus(:,1) == to);
                C_t(i, from_i) = 1;
                C_f(i, to_i) = 1;     
            end

            % make Pbusshift and Bbus
            obj.P_busshift = (C_f - C_t)' * obj.P_fshift;
            obj.B_f = diag(B_ff) * (C_f - C_t);
            obj.B_bus = (C_f - C_t)' * obj.B_f;
            
            % create Bbustildeinv
            B_bustilde = obj.B_bus(1:end-1, 1:end-1);
            obj.B_bustildeinv = inv(B_bustilde);
            
            % create switching matrix for generators C_g
            obj.C_G = zeros(obj.N_b, obj.N_G);
            for i = 1:obj.N_G
                bus_i = mpc.gen(i, 1);
                j = mpc.bus(:,1) == bus_i;
                obj.C_G(j, i) = 1;
            end
            
            % line ratings
            obj.P_fmax = mpc.branch(:, 6);   % Line ratings [MVA]
            
            % if all maximum line ratings are 0, make the line rating
            % unbounded (1000000000 MW)
            if sum(obj.P_fmax) == 0
                obj.P_fmax = 1e9 * ones(obj.N_l, 1);
            end

            % Create load curve based on values from the cases
            t = 1:24;
            P_Dmax = mpc.bus(:, 3);        % Demanded real power [MW]
            obj.P_D = P_Dmax;
%             obj.P_D = 0.5*P_Dmax*ones(1,24) + 0.5*P_Dmax*sin(t/28*pi) + ...
%                       P_Dmax*0.1*(rand(1,24)-0.5);
            
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
            obj.c_us = (obj.c_qu.*(obj.P_Gmax.^2)+obj.c_li.*obj.P_Gmax) ...
                                                        ./ obj.P_Gmax;
                                                    
            % downspinning costs are the less than upspinning costs
            obj.c_ds = 0.8*obj.c_us;

        end
        
        function set_WPG_bus(obj, WPG_bus_index)
            % set_WPG_bus(bus_index)
            % Sets the WPG bus to a new index, updates C_w
            %
            % Parameters
            % bus_index : the index of the bus infeed (scalar or vector)
            %
            % Returns
            % Nothing
            
            % determine number of buses
            obj.N_w = length(WPG_bus_index);
            
            % fill connection matrix for wind C_w
            obj.C_w = zeros(obj.N_b, obj.N_w);
            for i = 1:length(WPG_bus_index)
                obj.C_w(WPG_bus_index(i), i) = 1;
            end
            
        end
        
        % draws the network layout
        function draw_network(obj)
            
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
            
            % plot the graph and highlight special nodes
            G = graph(adj_matrix, labels{1});

            clf
            h = G.plot();
            highlight(h, gens, ...
                    'Marker', 's', 'MarkerSize', 7, 'NodeColor', 'r');
            highlight(h, winds, ...
                    'NodeColor', 'g', 'Marker', 'p', 'MarkerSize', 7);
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);

        end        
    end
end


