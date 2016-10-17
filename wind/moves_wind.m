function [P_twWP P_so] = moves_wind(Po,N_t,N_w,P_tr,P_s,N_chain)
% Po : initial state of wind power in p.u.
% N_t : number of states in the simulated chain (horizon)
% N_w : wind power scenarios 
% P_tr : Wind power transition matrix
% P_s : vector with elements the N_g states of the
% chain (in p.u. values)
% P_so : initial state of the chain (interlpolation in the grid)
% P_twWP : wind power scenario matrix (timesteps x scenarios)
% N_chain : number of elements of the chain

%% generate wind scenarios based on MCMC
% construct the cumulative probability transition matrix Pcum_tr
N_g = N_chain; % number of states of the chain, including '0' and 'P_N'
Pcum_tr = zeros(N_g,N_g);
Pcum_tr(:,1) = P_tr(:,1);
for i = 2:1:N_g
    Pcum_tr(:,i) = Pcum_tr(:,i-1) + P_tr(:,i);
end

P_twWP = []; % create the wind power scenario matrix
P_twWP_so = []; % monitor the initial states of the chain for each scenario
for n_sim = 1:1:N_w
    % map Po to a state of the chain
    [val_min i_so] = min(abs(P_s-Po)); % i_so is the index of the current state of the chain
    P_so = P_s(i_so); % value of the wind power at the current state of the chain
    chain = i_so; % states of the simulated chain
    Psim = []; % generated states for the wind power
    v = []; % store the selected random number u
    for i = 1:1:N_t
        u = rand;
        [val_min i_s] = min(abs(Pcum_tr(i_so,:)-u));
        if Pcum_tr(i_so,i_s)>=u
            chain = [chain i_s];
            Psim = [Psim P_s(i_s)];
            i_so = i_s;
        else
            chain = [chain i_s+1];
            Psim = [Psim P_s(i_s+1)];
            i_so = i_s + 1;
        end
        v = [v u];
    end
    P_twWP_so = [P_twWP_so P_so'];
    P_twWP = [P_twWP Psim'];
    clear Psim chain v
end

