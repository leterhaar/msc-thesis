function check_class(vars, types)
% check if incoming arguments are of the right type and dimensions, throws
% error if not the case

    % if it is not a cell, compare all arguments
    if not(isa(vars, 'cell'))
        assert(isa(vars, types), 'Not a %s', types);
    else
        assert(length(vars) == length(types));
        % loop over all arguments
        for i = 1:length(vars)
            assert(isa(vars{i}, types{i}), 'Not a %s', types{i});
        end
    end
end