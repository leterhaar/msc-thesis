function [Psim_case0f Psim_case0 Psim_case1] = moves_scenarios_forecast(N_t,Nstep,P_tr,P_s,N_w,P_b,i,N_chain,horizon)
% MOVES_SCENARIO_FORECAST generates scenarios of the wind power
% INPUT:
% N_t:      control horizon
% N_step:   optimization horizon
% P_tr:     transition matrix
% P_s:      vector with elements the N_g states of the chain
% N_W:      number of scenarios to generate
% P_b:      base value to scale the wind power
% i:        number of the day
% N_chain:  number of elements of the chain (needed in moves_wind.m)
% horizon:  form of wind prediction
% OUTPUT:
% Psim_case0f: forecasted wind power
% Psim_case0:  actual wind power
% Psim_case1:  generated scenarios 

plot_scenarios = 0; 
switch horizon
    case 1
        load WP_data_DA_24h
        Nlost = 16; % no action period for the day ahead market
    case 2
        load WP_data_MH_01h
        Nlost = 1;  % prediction horizon for moving horizon data
    case 3
        load WP_data_MH_04h
        Nlost = 4;  % prediction horizon for moving horizon data
    case 4
        load WP_data_MH_08h
        Nlost = 8;  % prediction horizon for moving horizon data
    otherwise
        disp('Not a valid case! In moves_scenarios_forecast.m');
end
load scale
%% generate wind power forecast + base (measured) case
% Matrix (hours x days) = (24 x 90)
day_h = 24;
offset_h = 1; % so that size(P_f)/day is an integer
Pf_new = Pf_new./scale;
Pm_new = Pm_new./scale;
for k = 1:1:N_t
    P_f_hour(k,:) = Pf_new(k:day_h:end-offset_h);
    P_m_hour(k,:) = Pm_new(k:day_h:end-offset_h);
    P_e_hour(k,:) = P_f_hour(k,:)-P_m_hour(k,:); 
end
% Defining real and forecast for day "i"
Psim_case0 = P_m_hour(:,i); % the "assumed" real scenario
Psim_case0f = P_f_hour(:,i); % the forecast

% Defining real and forecast for day "i-1"
Pdm1_case0 = P_m_hour(:,i-1); % the "assumed" real scenario
Pdm1_case0f = P_f_hour(:,i-1); % the forecast

switch horizon
    %% generate N_w scenarios for starting at t=1 and keep the parts t=Nlost+1:Nlost+Nstep (case 1)
    case 1 % Day Ahead Scenarios
        % Initial Forecast-Error at 24h-Nlost in day "i-1"
        Po = Pdm1_case0f(N_t-Nlost+1)-Pdm1_case0(N_t-Nlost+1);
        Psim_case1 = [];
        [Psim P_so] = moves_wind(Po,Nstep+Nlost,N_w,P_tr,P_s,N_chain);
        Psim = [Po*ones(1,N_w); Psim(1:end-1,1:end)];
        Psim_case1 = repmat(Psim_case0f,1,N_w)-Psim(Nlost+1:end,:);
        pos1 = find(Psim_case1(:)<0);
        Psim_case1(pos1) = 0;

        % Only to visualize (inclusive the part that will be cut off):
        if plot_scenarios
        Pall_case1 = repmat([Pdm1_case0f(end-Nlost+1:end);Psim_case0f],1,N_w)-Psim;
        pos1 = find(Pall_case1(:)<0);
        Pall_case1(pos1) = 0;
        end
        
    %% generate N_w scenarios for each timestep Nstep (cases 2,3,4)
    case {2,3,4}
        % Merge day "i" and day before "i-1"
        P2d_case0 = [Pdm1_case0;Psim_case0];
        P2d_case0f = [Pdm1_case0f;Psim_case0f];
        for j = N_t+1:1:2*N_t
            % Initial Forecast-Error:
            Po = P2d_case0f(j-Nlost,1)-P2d_case0(j-Nlost,1);
            [Psim P_so] = moves_wind(Po,Nlost+1,N_w,P_tr,P_s,N_chain);
            Psim = [Po*ones(1,N_w); Psim(1:end-1,1:end)];
            Psim_m = repmat(P2d_case0f(j-Nlost:j,1),1,N_w)-Psim;
            pos_m = find(Psim_m(:)<0);
            Psim_m(pos_m) = 0;
            Psim_case3_temp(j-N_t,:) = Psim_m(end,:);
            
            % Only to visualize (inclusive the part that will be cut off):
            if j-N_t == 1 && plot_scenarios
                Pall0_case3 = Psim_m;        
            elseif j-N_t == 16 && plot_scenarios
                Pall15_case3 = Psim_m;        
            end
        end
        Psim_case1 = Psim_case3_temp; 
    otherwise
        disp('Not a valid case! In moves_scenarios_forecast.m');
end
% convert from p.u -> MW
Psim_case0f = Psim_case0f*P_b;
Psim_case0 = Psim_case0*P_b;
Psim_case1 = (Psim_case1*P_b);

%% Plotting
if plot_scenarios
Pdm1_case0f = Pdm1_case0f*P_b;
Pdm1_case0 = (Pdm1_case0*P_b);
switch horizon
    case 1
        Pall_case1 = (Pall_case1*P_b);
        figure();title(['Forecast:',horizon,' - Day i=',num2str(i),' and Day i-1=',num2str(i-1),' (grey background)']); 
        hold on;grid on;box on;set(gca,'xtick',0:2:47,'xticklabel',[[0:2:23],[0:2:23]],'xlim',[0 47],'Layer','top');
        a = area([0 24], [10 10]);
        for i=1:1:size(Psim_case1,2)
            sc_co = plot((N_t-Nlost):1:2*N_t-1,Pall_case1(:,i),'m');     
            sc = plot((N_t):1:2*N_t-1,Psim_case1(:,i),'g');       
        end
        ac = plot(0:1:47,[Pdm1_case0;Psim_case0],'b');
        fc = plot(0:1:47,[Pdm1_case0f;Psim_case0f],'r');
        set(a,'YData',[max(get(gca,'YLim')),max(get(gca,'YLim'))],'Edgecolor','none','Facecolor',[0.9 0.9 0.9])
        legend([ac,fc,sc(1),sc_co(1)],'actual','forecast','scenarios','lost-scenarios');
    case {2,3,4}
        Pall0_case3 = (Pall0_case3*P_b);
        Pall15_case3 = (Pall15_case3*P_b);
        figure();title(['Forecast:',horizon,' - Day i=',num2str(i),' and Day i-1=',num2str(i-1),' (grey background)']); 
        hold on;grid on;box on;set(gca,'xtick',0:2:47,'xticklabel',[[0:2:23],[0:2:23]],'xlim',[0 47],'Layer','top');
        a = area([0 24], [10 10]);
        for i=1:1:size(Psim_case3,2)
            sc = plot(N_t:1:2*N_t-1,Psim_case3(:,i),'g');
            sc_co0 = plot(N_t-Nlost:1:N_t,Pall0_case3(:,i),'m');
            sc_co15 = plot(N_t+15-Nlost:1:N_t+15,Pall15_case3(:,i),'m');            
        end
        ac = plot(0:1:47,[Pdm1_case0;Psim_case0],'b');
        fc = plot(0:1:47,[Pdm1_case0f;Psim_case0f],'r');
        set(a,'YData',[max(get(gca,'YLim')),max(get(gca,'YLim'))],'Edgecolor','none','Facecolor',[0.9 0.9 0.9])
        legend([ac,fc,sc(1),sc_co0(1),sc_co15(1)],'actual','forecast','scenarios','lost-scenarios for h=0','lost-scenarios for h=15');
    otherwise
        disp('Not a valid case! In moves_scenarios_forecast.m');
end
end