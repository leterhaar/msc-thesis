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

clear agents
for ii = 1:m
    agents(ii) = agent(ac, wind, union(buses{ii,:}), buses{ii, ii}, 100);
end
p = progress('Running multi-area ADMM', 1000);
for it = 1:1000
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

%% plot residuals
res = 0;
for ii = 1:m
    res = res + agents(ii).residuals;
end

semilogy(res, 'linewidth', 2);
hold on; grid on; 
xlabel('iterations')
ylabel('residuals');
title('1000 iterations of 2 areas OPF');