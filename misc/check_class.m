function check_class(vars, types)
% check if incoming arguments are of the right type, throws
% error if not the case

    % if it is not a cell, compare all arguments
    if not(isa(vars, 'cell'))
        check_class_single(vars, types)
    else
        assert(length(vars) == length(types));
        % loop over all arguments
        for i = 1:length(vars)
            check_class_single(vars{i}, types{i})
        end
    end
end

% check type
function check_class_single(var, types)
    types_splitted = strsplit(types, '|');
    num_its_not = 0;
    for type = types_splitted
        if not(isa(var, type{:}))
            num_its_not = num_its_not + 1;
        end
    end
    if length(types_splitted) == num_its_not
        error('Not a %s', types);
    end
end

    