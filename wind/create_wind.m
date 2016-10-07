clear all
clc

% Load data for the test power system:
sys = getSystemData('threebus');

% Names of the wind prediction form:
% da24 -> 1, mh01 ->2, mh04 ->3, mh08 ->4
wind = 1;

%% Numerical data
Nsim = 1; % number of simulations
N_t = 24; % control horizon of each problem
Nstep = 24; % opimization horizon

N_B = sys.N_B; % buses n=1...N_B
N_G = sys.N_G; % generators i=1...N_G
N_L = sys.N_L; % loads j=1...N_L
N_l = sys.N_l; % lines nl=1...N_l
n_x = 3*N_G + 2*N_L + N_B; % Number of states

%% 
% percentage of wind power integration
WP_perc = 20;

% load profile
load_prof = [sys.load_prof_day];

% wind power base
load WP_data_DA_24h
P_b = (WP_perc/100) * max(load_prof) * getScaleWP(wind,Nsim);

% define number of scenarios
Eps_perc = 10;
beta = 1e-3; % confidence level
epsilon = Eps_perc/100; % violation parameter
N_w = ceil((2/epsilon)*(log(1/beta) + n_x)); % number of wind scenarios
               
% Generate wind scenarios for the scenario approach:
[P_tr, P_s, P_f, P_m, N_chain] = moves_TrM_forecast(); % Create the Transition matrix P_tr
Perror = null(P_tr'-eye(size(P_tr,1)))/sum(null(P_tr'-eye(size(P_tr,1))));
for i = 2:1:Nsim+1 
    [Psim_case0f{i-1}, Psim_case0{i-1}, Psim_case1{i-1}] = ...
        moves_scenarios_forecast(N_t,Nstep,P_tr,P_s,N_w,P_b,i,N_chain,wind);
end
% 
% % Generate quantiles for the benchmark approach:
% min_value = epsilon/2;
% max_value = 1-epsilon/2;
% index_chainMin = find(cumsum(Perror)<=min_value,1,'last');
% index_chainMax = find(cumsum(Perror)>=max_value,1,'first');
% for i = 2:1:Nsim+1
%     for k = 1:1:N_t
%         Pmin_quant1{i}(k,1) = Psim_case0f{i}(k,1) - P_b*P_s(index_chainMin);
%         Pmax_quant1{i}(k,1) = Psim_case0f{i}(k,1) - P_b*P_s(index_chainMax);
%     end
%     Pmin_quant1{i}(find(Pmin_quant1{i}<0)) = 0;
%     Pmax_quant1{i}(find(Pmax_quant1{i}<0)) = 0;
% end