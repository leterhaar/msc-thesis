classdef experiment
    
    properties
        
        name;
        network;
        wind;       
    end
    
    
    methods
        
        
        
        function exp = experiment(name, varargin)
            % prepares an experiment by creating a settings structure called
            % experiment and loads the wind and model networks
            %
            % PARAMETERS
            % ==========
            % network_name:           the name of the network (default case14)
            % network_formulation:    DC, P1, P2, P3 (default DC)
            % network_windbus:        bus id(s) of the wind bus of the network (default 9)
            % wind_N:               number of scenarios to use (default 100)
            % wind_Nm:              number of scenarios per agent (default 10)
            % wind_horizon:         number of timesteps (default 24)
            % wind_penetration:     percentage of wind penetration (default 0.2)
            % algorithm:            algorithm name (default ACC)

            %% construct settings variable
            exp.name = name;
            
            % reconstruct keys and values
            assert(rem(nargin-1,2) == 0, 'Please provide key-value pairs');

            fields = varargin(1:2:end);
            values = varargin(2:2:end);
            
            if ~isempty(fields) && strcmp(fields{1}, 'loading')
                return
            end
            
            
            % default experiment values
            exp.network = ...
                struct('name',              'case14',     ...
                       'formulation',       'AC',        ...
                       'windbus',           9);
            exp.wind = ...
                struct( 'N',                        100, ...
                        'horizon',                  24, ...
                        'penetration',              0.2, ...
                        'dummy',                    false);
            
            def_fields = fieldnames(exp);

            % loop over fields that need to be set
            for i = 1:length(fields)
                
                field_known = 0;
                if find(strcmp(def_fields, fields{i}))
                    exp.(fields{i}) = values{i};
                    field_known = 1;
                else

                    % see if it is a sublevel 
                    splitted = strsplit(fields{i}, {'.', '_'});
                    if length(splitted) == 2
                        first = splitted(1);
                        second = splitted(2);

                        % check if these fields exists
                        if find(strcmp(def_fields, first))
                            def_second_fields = fieldnames(exp.(first{:}));
                            if find(strcmp(def_second_fields, second))
                                exp.(first{:}).(second{:}) = values{i};
                                field_known = 1;
                            end
                        end
                    end
                end

                if not(field_known)
                    warning('Field "%s" is unknown. Typo?', fields{i});
                end
            end    
            
            % make sure that the paths are added
            if not(exist('AC_model', 'file'))
                addpath('../networks');
            end
            if not(exist('wind_model', 'file'))
                addpath('../wind');
            end
            

            if strcmpi(exp.network.formulation, 'DC')
                network = DC_model(exp.network.name);
                network.set_WPG_bus(exp.network.windbus);
                assignin('base', 'dc', network);
            else
                network = AC_model(exp.network.name);
                network.set_WPG_bus(exp.network.windbus);
                network.get_bags();
                assignin('base', 'ac', network);
            end
            exp.network.model = network;

            wind = wind_model(network, 24, exp.wind.penetration);
            if exp.wind.dummy
                wind.dummy(exp.wind.N);
            else
                wind.generate(exp.wind.N);
            end
            wind.shorter_horizon(exp.wind.horizon);
            exp.wind.model = wind;
            
        end
        
        function save(exp)
        %% saves the experiment
            save(['/Users/leterhaar/Dropbox/TU/Afstuderen/Ole/code/experiments/' exp.name], 'exp');
        end
        
        function disp(exp)
            disp(sprintf('''%s'': %s %s with %i scenarios over %i timesteps\n',...
                         exp.name, exp.network.formulation, exp.network.name, ...
                         exp.wind.N, exp.wind.horizon));
        end
        
        
    end
    
    methods(Static)
        
        function exp = load(name)
        %% load the experiment
            saved = load(['/Users/leterhaar/Dropbox/TU/Afstuderen/Ole/code/experiments/' name]);
            exp = experiment(name, 'loading', 1);
            
            % make sure that the paths are added
            if not(exist('AC_model', 'file'))
                addpath('../networks');
            end
            if not(exist('wind_model', 'file'))
                addpath('../wind');
            end
            
            % loop over properties and assign
            for propname = fieldnames(saved.exp)'
                exp.(propname{1}) = saved.exp.(propname{1});
            end
        end
    end
end