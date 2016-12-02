function [SDP] = read_sdpa(filename)
% reads an SDPA file and returns the system matrices

    % test if file exists
    if not(exist(filename, 'file') == 2)
        error('Can`t find %s anywhere...', filename);
    end

    % load file
    file = fopen(filename);

    % find the first line that does not start with a comment
    tline = fgetl(file);
    while tline(1) == '*' || tline(1) == '"'
        tline = fgetl(file);
    end

    % extract m
    m_str = strsplit(tline);
    m = str2num(m_str{1});

    % extract n 
    tline = fgetl(file);
    n_str = strsplit(tline);
    n = str2num(n_str{1});

    % get block sizes
    tline = fgetl(file);
    blocksizes_str = strsplit(tline, {' ', ',', '(', '{', ')', '}'});
    blocks = cell(m+1, 0);
    j = 1;
    N = 0;
    for i = 1:length(blocksizes_str)
        if not(isempty(blocksizes_str{i}))
            blocksize = abs(str2num(blocksizes_str{i}));
            [blocks{:, j}] = deal(zeros(blocksize));
            N = N + blocksize;
            j = j + 1;
        end
    end
    assert(n == j-1, 'block sizes not read correct');
    
    % get bs
    tline = fgetl(file);
    bs = arrayfun(@(x) str2double(x), strsplit(tline));
    SDP.bs = bs(not(isnan(bs)));
    
    assert(m == length(SDP.bs), 'no of constraints not matching');

    % fill blocks
    tline = fgetl(file);
    while ischar(tline)
       nline = arrayfun(@(x) str2double(x), strsplit(tline));
       matnum = nline(1)+1;
       blocknum = nline(2);
       i = nline(3);
       j = nline(4);
       entry = nline(5);

       blocks{matnum, blocknum}(i, j) = entry;
       if j ~= i
           blocks{matnum, blocknum}(j, i) = entry;
       end

       tline = fgetl(file);

    end

    As = cell(m+1,1);
    for i = 1:m+1
        As{i} = blkdiag(blocks{i, :});
    end

    % since in sdplib formulation, original problem is maximization
    SDP.C = -As{1};
    
    SDP.As = As(2:end);

    SDP.X = sdpvar(N);
end
