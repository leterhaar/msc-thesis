function [ScaleWP] = getScaleWP(horizon,Nsim)
% GETSCALEWP calculates a scale factor depending on mean(WP_act)
% INPUT: horizon: form of WP prediction
% Nsim: number of simulations
% OUTPUT: ScaleWP = scale factor
switch horizon
    case 1
        load WP_data_DA_24h
        ScaleWP = 23900/mean(Pm_new(25:(Nsim+1)*24));
    case 2
        load WP_data_MH_01h
        ScaleWP = 23900/mean(Pm_new(25:(Nsim+1)*24));    
    case 3
        load WP_data_MH_04h
        ScaleWP = 23900/mean(Pm_new(25:(Nsim+1)*24));
    case 4
        load WP_data_MH_08h
        ScaleWP = 23900/mean(Pm_new(25:(Nsim+1)*24));
    otherwise
        disp('Not a valid case! In getScaleWP.m');
end