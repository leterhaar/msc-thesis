function verify(that, varargin)
% debug function
    if not(that)
        if nargin < 2
            fprintf('\nVerification failed\n');
        else
            msg = sprintf(varargin{1}, varargin{2:end});
            fprintf('\nVerification failed: %s\n', msg);
        end
        dbstack
        evalin('caller', 'hasVerificationErrors = 1;');
    end
end