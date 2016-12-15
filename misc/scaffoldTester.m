function scaffoldTester(name)
% creates a scaffold for the test with the right name in place

    assert(not(isempty(name)), 'Can not be run with empty name');
    
    filename = sprintf('test%s%s.m', upper(name(1)), name(2:end));
    template = fileread('templateTests.m');
    
    
    fh = fopen(filename, 'w');
    fprintf(fh, template(:), filename(1:end-2));
    fclose(fh);
    
    fprintf('\nCreated %s\n', filename);
    edit(filename)
    
end
