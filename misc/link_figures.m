%% link_figures
% function that walks over all dirs makes symlinks 
%  to the .eps files in there

function link_figures(basedir, epsdir)

    if nargin < 1 
        basedir = '~/Dropbox/TU/Afstuderen/Ole/code';
    end
    
    if nargin < 2
        epsdir = '~/Dropbox/TU/Afstuderen/Ole/tex/main/report/img';
    end
    
    % call nested function
    loop_over_dir(basedir, epsdir);

end

function loop_over_dir(dirname, epsdir)
% loop_over_dir(dirname) 
    epssfound = 0;
    % read dir contents
    contents = dir(dirname);

    for file = contents'
        if logical(file.isdir) && isempty(strfind(file.name, '.'))
            loop_over_dir([dirname '/' file.name], epsdir)
        end
        
        % check if it is a .eps
        if not(isempty(strfind(file.name, '.eps')))
            
            target = [epsdir '/' file.name];
            
            % check if file already linked
            if not(exist(target, 'file'))
                system(['ln -s ' dirname '/' file.name ' ' target]);
                epssfound = epssfound + 1;
            end
        end
    end
    
    if epssfound > 0
        foldernames = strsplit(dirname, '/');
        last_two = strjoin(foldernames(end-1:end), '/');
        fprintf('Linked %i new file(s) in %s\n', epssfound, last_two)
    end
end
