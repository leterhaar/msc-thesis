function mv(varargin)
% alias for movefile
    movefile(varargin{:});

    openfiles = matlab.desktop.editor.getAll;

    % close the first file and open the second
    for i = 1:length(openfiles)
        parts = strsplit(openfiles(i).Filename, '/');
        if strcmp(parts(end), varargin{1})
            openfiles(i).close();
        end
    end
    
    edit(varargin{2});
end