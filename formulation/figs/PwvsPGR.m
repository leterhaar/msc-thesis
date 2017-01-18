addpath('../../misc');

%% extract Pgopt from Wsscen
for i = 1:N
    for j = 1:ac.N_G
        k = ac.Gens(j);
        P_Gopt(i,j) = trace(ac.Y_k(k)*W_sscen(:,:,i)) - ac.C_w(k)*wind.P_w(t, i) + ac.P_D(t, k);
    end
end

% plot
% P_Gstacked = cumsum(P_Gopt, 2);
initfig('P^G + R vs P^m');
plot(repmat(wind.P_m(t,:)', 1, ac.N_G), P_Gopt);
xlim([min(wind.P_m(t,:)) max(wind.P_m(t,:))]);
legend('G1','G2','G3','G4','G5')