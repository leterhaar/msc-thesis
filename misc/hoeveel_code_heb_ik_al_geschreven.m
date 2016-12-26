function hoeveel_code_heb_ik_al_geschreven(dirname)

    if nargin < 1 
        dirname = '~/Dropbox/TU/Afstuderen/Ole/code';
    end
    
    % call nested function
    bytes = loop_over_dir(dirname, 0);
    
    fprintf('\nTOTAL: %s\n', format_bytes(bytes));

end

function bytes = loop_over_dir(dirname, bytes_start)
% loop_over_dir(dirname) 
    % read dir contents
    bytes = bytes_start;
    files_found_in_this_folder = 0;
    contents = dir(dirname);

    for file = contents'
        if logical(file.isdir) && isempty(strfind(file.name, '.'))
            bytes = loop_over_dir([dirname '/' file.name], bytes);
        end
        
        % check if extension ends with .m
        if length(file.name) > 2
            if strcmp(file.name(end-1:end), '.m')
                files_found_in_this_folder = files_found_in_this_folder + 1;
                bytes = bytes + file.bytes;
            end
        end
    end
    if files_found_in_this_folder > 0 
        foldernames = strsplit(dirname, '/');
        last_two = strjoin(foldernames(end-1:end), '/');
        fprintf('%-25s %3i .m files: %s\n', ...
                            last_two, files_found_in_this_folder,...
                            format_bytes(bytes-bytes_start));
    end
end

function res = format_bytes(bytes)
    Mb = 2^20;
    Kb = 2^10;
    format = '%.2f';
    Mbs = floor(bytes / Mb);
    Kbs = floor(bytes / Kb);

    if Mbs > 0
        res = sprintf([format ' Mb'], bytes/Mb);
    elseif Kbs > 0
        res = sprintf([format ' Kb'], bytes/Kb);
    else
        res = sprintf([format ' b'], bytes);
    end
    
end




