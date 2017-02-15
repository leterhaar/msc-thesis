function init_experiment(varargin)
% prepares an experiment by creating a settings structure called
% experiment and loads the wind and network models
%
% PARAMETERS
% ==========
% model_name:           the name of the model (default case14)
% model_formulation:    DC, P1, P2, P3 (default DC)
% model_windbus:        bus id(s) of the wind bus of the model (default 9)
% wind_N:               number of scenarios to use (default 100)
% wind_Nm:              number of scenarios per agent (default 10)
% wind_horizon:         number of timesteps (default 24)
% wind_penetration:     percentage of wind penetration (default 0.2)
% algorithm:            algorithm name (default ACC)

    %% construct settings variable

    % reconstruct keys and values
    assert(rem(nargin,2) == 0, 'Please provide key-value pairs');

    fields = varargin(1:2:end);
    values = varargin(2:2:end);

    % default experiment values
    exp = struct( ...
        'model',    ...
        struct('name',              'case14',     ...
               'formulation',       'DC',        ...
               'windbus',           9),          ...
        'wind',             ...
        struct( 'N',                        100, ...
                'Nm',                       10, ...
                'horizon',                  24, ...
                'penetration',              0.2, ...
                'dummy',                    false), ...
        'algorithm',                'ACC');

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
    
    %% Create models
    addpath('../networks');
    addpath('../wind');
    addpath('../misc');
    addpath('../formulation');
    
    if strcmpi(exp.model.formulation, 'DC')
        model = DC_model(exp.model.name);
        model.set_WPG_bus(exp.model.windbus);
        assignin('base', 'dc', model);
    else
        model = AC_model(exp.model.name);
        model.set_WPG_bus(exp.model.windbus);
        assignin('base', 'ac', model);
    end
    exp.model.network = model;
    
    wind = wind_model(model, exp.wind.horizon, exp.wind.penetration);
    if exp.wind.dummy
        wind.dummy(exp.wind.N);
    else
        wind.generate(exp.wind.N);
    end
    assignin('base', 'wind', wind);
    exp.wind.model = wind;
    
    % decide number of agents and cut indices
    cut_indices = 1:exp.wind.Nm:exp.wind.N;
    m = length(cut_indices);
    assignin('base', 'm', m);
    
    % build cut indices: agent i goes from cut(i,1) to cut(i,2)
    cut = zeros(m,2);
    for i = 1:m
        cut(i,1) = cut_indices(i);
        if i < m
            cut(i,2) = cut_indices(i+1)-1;
        else
            cut(i,2) = exp.wind.N;
        end
    end
    assignin('base', 'cut', cut);
    
    assignin('base', 'experiment', exp);

end