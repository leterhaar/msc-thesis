function check_class(vars, types)
% check if incoming arguments are of the right type, throws
% error if not the case

    % if it is not a cell, compare all arguments
    if not(isa(vars, 'cell'))
        assert(isa(vars, types), 'Not a %s', types);
    else
        assert(length(vars) == length(types));
        % loop over all arguments
        for i = 1:length(vars)
            types_splitted = strsplit(types{i}, '|');
            nums_its_not = 0;
            for type = types_splitted
                try
                    assert(isa(vars{i}, type{:}), 'Not a %s', types{i});
                catch
                    nums_its_not = nums_its_not + 1;
                end
            end
            if nums_its_not == length(types_splitted)
                assert(false, 'Not a %s', types{i});
            end
                
        end
    end
end