%% Printer class
% Prints status of loops and an estimate of the duration
% to the command line

classdef progress < handle

    properties
        
        N;              % the number of loops
        i;              % current index
        start_time;     % time of starting the loop
        previous_perc   % previous percentage
        previous_update % previous update time
    end

    methods
       
        function obj = progress(msg, N)
            % constructor function
            % prints the message to the command line
            % and starts the timer
            obj.N = N;
            obj.i = 0;
            obj.start_time = tic;
            obj.previous_perc = -1;
            obj.previous_update = tic;

            fprintf('%s: ', msg);
            obj.print_percentage(0);
        end

        function print_percentage(obj, perc)
            % prints the percentage bar
            dec = floor(perc/5);
            
            % only update after a second to prevent flickering
            if (toc(obj.previous_update) > 0.5 && obj.previous_perc ~= floor(perc)) || obj.i == obj.N || obj.i == 0
                obj.previous_update = tic;
                obj.previous_perc = floor(perc);
                % delete the old bar
                if obj.i ~= 0
                    fprintf(repmat('\b', 1, 52));
                end
                
                % running bar
                if obj.i < obj.N
                    arrow = [repmat('=',1,dec) '>' repmat('.',1,19-dec)];
                    fprintf('[%s] %3d%%   | ', arrow, floor(perc));

                    % print time
                    time_remaining = ceil(toc(obj.start_time) / ...
                                                  obj.i * (obj.N - obj.i));
                    if isinf(time_remaining)       
                        fprintf('??:??:?? remaining\t\n');
                    else
                        hours_remaining = floor(time_remaining / 3600);
                        time_remaining = rem(time_remaining, 3600);
                        minutes_remaining = floor(time_remaining / 60);
                        seconds_remaining = rem(time_remaining, 60);
                        fprintf('%02i:%02i:%02i remaining\t\n', hours_remaining, minutes_remaining, seconds_remaining);
                    end
                
                % complete
                else
                    fprintf('done in %.3f seconds\r', toc(obj.start_time));
                end
            end
        end
        
        function ping(obj)
            % updates the printer
            obj.i = obj.i + 1;
            perc = min(obj.i / obj.N * 100, 100);
            obj.print_percentage(perc);
        end
    end
end
