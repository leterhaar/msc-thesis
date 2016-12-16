function verify(that, msg, dostop)
% debug function
    if not(that)
        if nargin < 2
            fprintf('\nVerification failed\n');
        else
            fprintf('\nVerification failed: %s\n', msg);
        end
        dbstack
        if nargin > 2
            evalin('caller', 'keyboard');
        end
        evalin('caller', 'hasVerificationErrors = 1;');
    end
end