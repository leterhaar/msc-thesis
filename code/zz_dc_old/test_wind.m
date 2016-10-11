addpath('../Wind Model/')
m = DC_model('case30');
w = wind_model(m, 24);

epsilon = 1e-1;
zeta = m.N_G;
beta = 1e-5;

N_w = ceil(2/epsilon*(zeta-1+log(1/beta)));

for i = 1:20
    [P_wf, P_w, P_wscen] = w.simulate(N_w);
    clf
    t = 1:24;
    grid on
    hold on
    plot(t, P_wscen, ':g')
    plot(t, P_wf, 'LineWidth', 3)
    plot(t, P_w, 'LineWidth', 3)
    pause
end
