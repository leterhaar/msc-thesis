function test(name)
% init test variables
    yalmip('clear');
    warning on
    % create testCase object for verification
    close all
    for file = dir'
        if strcmp(file.name,'test_init.m')
            path = strsplit(pwd, '/');
            fprintf('Initializing tests in %s\n', path{end});
            eval('test_init');
            fprintf('   OK\n');

            if nargin == 1
                test_sequence = {name};
            end

            for i = randperm(length(test_sequence))
                all_ok = 1;
                scriptname = sprintf('test_%s', test_sequence{i});
                fprintf('Running %s\n   ', scriptname)
                try
                    hasVerificationErrors = 0;
                    eval(scriptname);
                    if not(hasVerificationErrors)
                        fprintf('OK\n');
                    else
                        fprintf('Failed\n');
                        all_ok = 0;
                    end
                catch ME
                    fprintf('Failed \n\n%s\n', ME.getReport());
                    dbstack;
                    keyboard;
                    all_ok = 0;
                end

            end
            
            yalmip('clear');
            
            if all_ok
                fprintf('All test OK\n');
            else
                fprintf('Some test failed\n');
            end
        end
    end
end