% Main function to generate tests
function tests = testVacuum
    tests = functiontests(localfunctions);
end

% Setup function
function setupOnce(testCase) 
    % make subdir and enter
    testdir = randomString(10);
    testCase.TestData.testdir = testdir;
    mkdir(testdir);
    copyfile('vacuum.m', [testdir '/vacuum.m']);
    cd(testdir);
    
    % create 5 random empty files and non empty files
    for i = 1:5
        empty_filename = randomString(10);
        testCase.TestData.empty{i} = empty_filename;
        fh = fopen(empty_filename, 'w');
        % .. empty file....
        fclose(fh);
        
        not_empty_filename = randomString(10);
        fh = fopen(not_empty_filename, 'w');
        testCase.TestData.not_empty{i} = not_empty_filename;
        % add some random text
        fprintf(fh, '%s', randomString(20));
        fclose(fh);
    end
end

function testFunction(testCase)
    deleted_files = vacuum;
    
    verifyEqual(testCase, length(deleted_files), 5);
    search_pattern = ['(' strjoin(testCase.TestData.empty, '|') ')'];
    for i = 1:length(deleted_files)
        verifyMatches(testCase, deleted_files{i}, search_pattern);
    end
    
end

function teardownOnce(testCase)
    cd ..
    rmdir(testCase.TestData.testdir, 's');
end

function string = randomString(N)
    alphabet = char(['a':'z' 'A':'Z' '0':'9']);
    n = length(alphabet);
    string = alphabet(randperm(n, N));
end