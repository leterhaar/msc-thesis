classdef subnetwork < AC_model
    
    properties
        
        buses;
        
    end
    
    
    methods
        
        function obj = subnetwork(ac, all_buses, own_buses)    
        % sac = subnetwork(ac, buses)    
        %
        % creates an AC_model with all data matrices that only have entries
        % that correspond to the buses
            
            % store bus set
            obj.buses = all_buses;
            
            % find lines in subnetwork
            lines = [];
            obj.from_to = [];
            for k = 1:ac.N_l
                from = ac.from_to(k, 1);
                to = ac.from_to(k, 2);
                if ismember(from, all_buses) && ismember(to, all_buses)
                    new_from = find(from == all_buses);
                    new_to = find(to == all_buses);
                    obj.from_to = [obj.from_to; new_from new_to];
                    lines = [lines; k];
                end
            end
            
            
            % find gens in subnetwork
            obj.C_G = ac.C_G(own_buses, :);
            obj.Gens = find(sum(obj.C_G) > 0);
            obj.C_G = obj.C_G(:, obj.Gens);
            
             % loop over properties and assign small values
            for propname = fieldnames(ac)'
                key = propname{1};
                value = ac.(key);
                
                % only select the buses
                if ismember(key, {'P_max','P_min','Q_max','Q_min','V_min','V_max', 'C_w'}) 
                    obj.(key) = value(all_buses, :);

                % admittance matrices
                elseif ismember(key, {'Y','Ystar'})
                    obj.(key) = value(all_buses, all_buses);
                
                % line limits
                elseif ismember(key, {'S_max', 'P_lmmax', 'y_sh'})
                    obj.(key) = value(lines, :);
                
                % load profiles
                elseif ismember(key, {'P_D', 'Q_D'})
                    obj.(key) = value(:, all_buses);
                
                % cost vectors
                elseif ismember(key, {'c_qu', 'c_li', 'c_us', 'c_ds'})
                    obj.(key) = value(obj.Gens);
                    
                else
                    %fprintf('%s not stored \n',key);
                end
            end
            
            if not(ismember(ac.refbus, all_buses))
                obj.refbus = [];
            else
                obj.refbus = ac.refbus;
            end            

            obj.C_w = ac.C_w(all_buses, :);
            obj.N_w = sum(sum(ac.C_w));
            obj.N_b = length(all_buses);
            obj.N_G = length(obj.Gens);
            obj.N_l = size(obj.from_to, 1);
            
            obj.model_name = sprintf('%s-sub', ac.model_name);
            
            % set P_max, P_min and Q_max, Q_min to NaN for any buses outside own network
            neighbours = not(ismember(all_buses, own_buses));
            obj.P_max(neighbours) = NaN;
            obj.P_min(neighbours) = NaN;
            obj.Q_max(neighbours) = NaN;
            obj.Q_min(neighbours) = NaN;
            
        end
        
        function disp(obj)
            fprintf('An sub-network from %s, only using buses %s\n\n', obj.model_name(1:end-4), num2str(obj.buses));
        end
            
    end 
    
    
end
    