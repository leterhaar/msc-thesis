%% Constructor function
% loads transition matrix and prepares wind simulator
% 
% Arguments
% ---------
% m : network model (to extract load profile)
% N_t : number of time steps
% wind_penetration : the share of power which is wind 
%
% Returns:
% --------
% nothing

classdef wind_model < handle
    
    properties
        N_t;        % control horizon
        N_step;     % optimization horizon
        P_tr;       % transition matrix
        P_s;        % vector with elements the N_g states of the chain
        P_b;        % base value to scale the wind power
        i;          % number of the day
        N_chain;    % number of elements of the chain (needed in moves_wind.m)
        horizon;    % form of wind prediction
        P_w;        % wind power scenarios
        P_wf;       % forecasted wind power
        P_m;        % wind power mismatch
    end
   
    
    methods
        
        function obj = wind_model(m, N_t, wind_penetration)
        %% Constructor function
        % loads transition matrix and prepares wind simulator
        % 
        % Arguments
        % ---------
        % m : DC_model (to extract load profile)
        % N_t : number of time steps
        % wind_penetration : the share of power which is wind 
        %
        % Returns:
        % --------
        % nothing

            if nargin > 0
                % share of WPG in network
                if nargin < 3
                    wind_penetration = 0.2;
                end
                % wind power base
                obj.P_b = wind_penetration * max(max(m.P_D)) * getScaleWP(1,1);

                % define transition matrix
                [obj.P_tr, obj.P_s, ~, ~, obj.N_chain] = moves_TrM_forecast();

                obj.N_t = N_t;
                obj.N_step = N_t;
                obj.i = 2;
                obj.horizon = 1;
            end
        end
        
        function [P_wf, P_w] = generate(obj, N_w, store)
            %% [P_wf, P_w, P_wscens] = generate(obj, N_w)
            % simulate the wind and return scenarios
            %
            % Parameters: 
            % N_w : number of scenarios to generate
            % store: bool, set true for storing in class. default 1
            %
            % Returns:
            % P_wf : forecasted wind power (N_t x 1)
            % P_w : wind trajectory scenarios (N_t x N_w)

            if nargin < 3
                store = 1;
            end
            
            [P_wf, ~, P_w] = ...
                moves_scenarios_forecast(obj.N_t, obj.N_step, obj.P_tr, ...
                                         obj.P_s, N_w, obj.P_b, ...
                                          obj.i, obj.N_chain, 1);    
                                      
            if store
                obj.P_wf = P_wf;
                obj.P_w = P_w;
                obj.P_m = P_w - (P_wf * ones(1, N_w));
            end
        end
        
        function plot(obj, t)
        % PLOT plots the scenarios
        
            dock
            clf
            hold on
            grid on
            
%             h = area(1:24, [min(obj.P_wscens, [], 2) ...
%                             max(obj.P_wscens, [], 2)]);
%             h(1).FaceColor = 'none';
%             h(1).EdgeColor = 'none';
%             h(2).EdgeColor = 'none';
%             h(2).FaceColor = ones(3,1)*0.95;
%           
            h = plot(-1, -1, 'LineWidth', 3);
            blue = get(h, 'Color');
            plot(obj.P_w, '--g');
            plot(obj.P_wf, '-', 'LineWidth', 3, 'Color', blue);
            
            % plot the time
            if nargin > 1
                
                plot(t, obj.P_wf(t), 's', 'markerfacecolor', blue, ...
                                          'markeredgecolor', blue);
                                      
                % loop over scenarios
                for Pw = obj.P_w
                    
                    % change arrow direction to illustrate up- and
                    % downspinning
                    if Pw(t) > obj.P_wf(t)
                        plot(t, Pw(t), 'v', 'markerfacecolor', 'g', ...
                                            'markeredgecolor', 'g');
                    else
                        plot(t, Pw(t), '^', 'markerfacecolor', 'g', ...
                                            'markeredgecolor', 'g');                        
                    end
                end

                
            end
            
            xlim([1 24]);
            xlabel('Time [h]');
            ylabel('Wind power [p.u.]');
            title('Wind power forecast and scenarios');
            legend('Forecasted', 'Scenarios');
        end
        
        function use_forecast(obj)
        %% use only the forecast value
            assert(~isempty(obj.P_w), 'Please generate scenarios first');
            [~,N_w] = size(obj.P_w);
            obj.P_w = obj.P_wf * ones(1, N_w);
            obj.P_m = obj.P_m.*0;   
            
        end
        
        function use_extremes(obj, t)
        %% use only the minimal and maximal forecast value
            assert(~isempty(obj.P_w), 'Please generate scenarios first');
            
            [~, min_id] = min(obj.P_w(t, :));
            [~, max_id] = max(obj.P_w(t, :));
            obj.P_w = obj.P_w(:, [min_id, max_id]);
            obj.P_m = obj.P_m(:, [min_id, max_id]);
        end
        
        function dummy(obj, N, min, max)
        %% creates a dummy wind variable to compare results
            if nargin < 4
                max = 2;
            end
            if nargin < 3
                min = 0;
            end
            
            obj.P_wf = ones(obj.N_t, 1) * obj.P_b * 0.2;
            obj.P_w = obj.P_wf * linspace(min, max, N);
            obj.P_m = obj.P_w - (obj.P_wf * ones(1, N));
        end
        
        function s = slice(obj, i_start, i_end)
        %% returns a slice of the generated wind
            if nargin < 3
                tmp = i_start;
                i_start = tmp(1);
                i_end = tmp(2);
            end
            s = struct( 'P_wf', obj.P_wf, ...
                        'P_w', obj.P_w(:, i_start:i_end), ...
                        'P_m', obj.P_m(:, i_start:i_end));
        end
    end
 
end