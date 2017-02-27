clear
yalmip('clear');
t = 1;
ac = AC_model('case14');
wind = wind_model(ac, 24, 0.2);
wind.generate(1);

bus1 = [1 2 3 4 5 7 8];
bus2 = [6 9 10 11 12 13 14];
overlap = [4 5 6 7 9];
buses = {bus1, overlap; overlap, bus2};
m = 2;
mus = [10 20 50 100 500];
N = 100;
res = zeros(5, N);
for mu_it = 1:5
    mu = mus(mu_it);
    clear agents
    for ii = 1:m
        agents(ii) = agent(ac, wind, union(buses{ii,:}), buses{ii, ii}, mu);
    end
    p = progress('Running multi-area ADMM', N);
    for it = 1:N
        for ii = 1:m
            agents(ii).update();
        end

        for ii = 1:m
            for jj = 1:m
                if ii == jj
                    continue
                end
                the_buses = buses{ii, jj};
                agents(ii).receive(agents(jj).broadcast(the_buses), the_buses);
            end
        end
        p.ping();
    end
    % aggregate residuals
    for ii = 1:m
        res(mu_it, :) = res(mu_it, :) + agents(ii).residuals;
    end
end
%%
clf
semilogy(res', 'linewidth', 2);
grid on
xlabel('iterations')
ylabel('residuals');
title('Effect of \mu on multi area ADMM convergence');
legend(cellfun(@(x) ['\mu = ' x], strsplit(num2str(mus)), 'UniformOutput', 0));

%% compare with centralized solution

Wf = sdpvar(2*ac.N_b);
Cons = [feasibleW(Wf, wind.P_wf, ac); Wf >= 0];
Obj = objective_PG(Wf, ac, wind);

optimize(Cons, Obj);
Wc = value(Wf);
[~, Vc] = completeW(Wc);

X = [real(Vc); imag(Vc)];
Wc = X * X';