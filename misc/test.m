function test(name)
% init test variables

if exist('test_init.m', 'file') == 2
    fprintf('Initializing tests\n');
    eval('test_init');
    fprintf('   OK\n');
end

if nargin == 1
    test_sequence = {name};
end

for i = 1:length(test_sequence)
    all_ok = 1;
    scriptname = sprintf('test_%s', test_sequence{i});
    fprintf('Running test_%s\n   ', scriptname)
    try
        eval(scriptname);
        fprintf('OK\n');
    catch ME
        fprintf('\nFailed: \n%s\n', ME.getReport());
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