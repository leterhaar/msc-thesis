classdef agent < handle
     
    properties
       % problem data
       ac; wind;
       buses; mu;
       
       % selection of entries from neighouring submatrices
       s_n; 
       
       % auxiliary variables and multipliers
       tildeW_n; Lambda_n;
       
       % local matrix
       W;
       
       % residuals
       residuals;
    end
    
       
    methods 
        
        function obj = agent(ac, wind, all_buses, own_buses, mu)
        % constructs an agent with its own sub-model and sub-wind
        % realization
        
            obj.ac = subnetwork(ac, all_buses, own_buses);
            obj.wind = wind;
            obj.buses = all_buses;
%             obj.own_buses = own_buses;
            obj.mu = mu;
            obj.W = zeros(2*obj.ac.N_b);
            obj.W(1:obj.ac.N_b, 1:obj.ac.N_b) = ones(obj.ac.N_b);
            obj.Lambda_n = cell(0,1);
            obj.residuals = [];
        end
        
        function update(obj)
            
            yalmip('clear');
            
            % evaluate the timestep in the caller 
            % NB this is a hack to work with one hour problems
            t = evalin('caller', 't');            
            
            % formulate optimization problem
            Constraints = [];
            Objective = [];
            Wf = sdpvar(2*obj.ac.N_b);
            
            % feasible W and PSD
            Constraints = [feasibleW(Wf, obj.wind.P_wf, obj.ac), Wf >= 0];
            Objective = objective_PG(Wf, obj.ac, obj.wind);
            
            for ii = 1:length(obj.s_n)
                Objective = Objective + obj.mu / 2 * norm( ...
                      Wf(obj.s_n{ii}, obj.s_n{ii}) - ...
                      obj.tildeW_n{ii} + ...
                      obj.Lambda_n{ii} / obj.mu, 'fro');
            end
           
            % optimize
            status = optimize(Constraints, Objective, sdpsettings('verbose', 0));
            verify(not(status.problem), status.info);

            % store solution
            obj.W = value(Wf);
            
            % store residuals
            obj.residuals(end+1) = 0;
            for ii = 1:length(obj.s_n)
                obj.residuals(end) = obj.residuals(end) + ...
                                     norm(obj.W(obj.s_n{ii}, obj.s_n{ii}) - ...
                                          obj.tildeW_n{ii}, 'fro');
            end
            
            % reset all auxiliary variables and submatrices
            obj.s_n = cell(0, 1);
            obj.tildeW_n = cell(0, 1);            
        end
        
        function W_ba = broadcast(obj, buses_neighbour)
            % returns the part of the local matrix that corresponds to
            % buses neighbour
            
            % find which local indices correspond to theses buses
            n = length(buses_neighbour);
            select = nan(2*n,1);
            for ii = 1:n
                select(ii) = find(buses_neighbour(ii) == obj.buses);
            end
            select(n+1:end) = select(1:n)+obj.ac.N_b;
            
            % broadcast those
            W_ba = obj.W(select, select);
        end
        
        function receive(obj, W_ab, buses_received)
            % stores the received submatrices and the corresponding entries
            % in local matrices, and stores the auxiliary variables and
            % multipliers corresponding to the submatrix
            
            % find which buses these correspond to in local matrix
            n = length(buses_received);
            s = nan(2*n,1);
            for ii = 1:n
                s(ii) = find(buses_received(ii) == obj.buses);
            end
            s(n+1:end) = s(1:n)+obj.ac.N_b;
            
            % store W_n and selection vector s_n
            obj.s_n{end+1} = s;
            
            % update auxiliary variables
            obj.tildeW_n{end+1} = 0.5*(obj.W(s, s) + W_ab);
            
            % update multipliers
            k = length(obj.s_n);
            if length(obj.Lambda_n) < k
                obj.Lambda_n{k} = zeros(length(s)); % initiate with zeros
            end
            
            obj.Lambda_n{k} = obj.Lambda_n{k} + ...
                              obj.mu / 2 * (obj.W(s, s) - W_ab);
            
        end
        
    end
end
