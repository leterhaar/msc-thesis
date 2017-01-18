function ls
while true
    clc
    display(pwd)
    display('---------------------------------------');
    list = dir;
    dirs = find([list.isdir] == 1);
    display('|-> [0] (mkdir)');

    for i = 1:length(dirs)-1
        display(['|-> [' int2str(i) '] ' list(dirs(i+1)).name]);
    end

    mfiles = getfield(what, 'm');
    [N,~] = size(mfiles);
    if(N > 0)
        display('|------ MFILES -------|');
    end
    for j = 1:N
        display(['|---> [' char('a'+j-1) '] ' mfiles{j}]);
    end
    models = getfield(what, 'slx');
    [Nm,~] = size(models);
    if(Nm > 0)
        display('|------ MODELS -------|');
    end
    for i = 1:Nm
        display(['|---> [' char('a'+j+i-1) '] ' models{i}]);
    end
    cmd = input('','s');
    if isempty(cmd)
        clc
        break
    end
    k = str2num(['int32(' cmd ')']);
    if ~isempty(k) && imag(k) == 0
        if k == 0
            clc
            newname = input('New folder name? ', 's');
            if ~isempty(newname)
                mkdir(newname);
                cd(newname);
            end
        elseif k+1 <= length(dirs);
            cd([pwd '/' list(dirs(k+1)).name '/']);
        end
    else
        k = cmd-96;
        if k <= N
            edit(mfiles{k})
            clc
            break;
        end
        if k-N <= Nm
            open(models{k-N})
            clc
            break;
        end
    end
end
end