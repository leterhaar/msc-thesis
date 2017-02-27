clear
yalmip('clear');
t = 1;
ac = AC_model('case14');
wind = wind_model(ac, 24, 0.2);
wind.generate(1);

bus1 = [1 2 3 4 5 7 8];
bus2 = [6 9 10 11 12 13 14];
overlap1 = [4 5 6 7 9];
overlap2 = [4 5 6 7 9];
buses = {bus1, overlap1; overlap2, bus2};
m = 2;
mu = 1000;
N = 250;
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

%% aggregate residuals
res = zeros(1,N);
for ii = 1:m
    res = res + agents(ii).residuals;
end
%%
clf
semilogy(res', 'linewidth', 2);
grid on
xlabel('iterations')
ylabel('residuals');
title('Multi area ADMM convergence');

%% compare with centralized solution

Wc = cell(m,1);
Cons = [];
Obj = [];

for ii = 1:m
    s_ac = subnetwork(ac, union(buses{ii,:}), buses{ii,ii});
    
    Wc{ii} = sdpvar(2*s_ac.N_b);
    
    Cons = [Cons, feasibleW(Wc{ii}, wind.P_wf, s_ac), Wc{ii} >= 0];
    Obj = Obj + objective_PG(Wc{ii}, s_ac, wind);
end

% overlapping constraints

for ii = 1:m
    for jj = ii:m
        if ii == jj
            continue
        end
        indices1 = [];
        indices2 = [];
        Nb1 = length(Wc{ii})/2;
        Nb2 = length(Wc{jj})/2;
        for bus = buses{ii, jj}
            index1 = find(bus == union(buses{ii, :}));
            index2 = find(bus == union(buses{jj, :}));
            indices1 = [indices1; index1; index1+Nb1];
            indices2 = [indices2; index2; index2+Nb2];
        end
        Cons = [Cons, Wc{ii}(indices1, indices1) == ...
                      Wc{jj}(indices2, indices2)];
    end
end

optimize(Cons, Obj)
Wc = values_cell(Wc);
min(check(Cons))
%% 
Vd = nan(ac.N_b, m);
for ii = 1:m
    
    % get full local matrix 
    Wd = Wc{ii};
    s_ac = subnetwork(ac, union(buses{ii, :}), buses{ii, ii});
    W = sdpvar(2*s_ac.N_b);
    Cons_d = [feasibleW(W, wind.P_wf, s_ac), W >= 0];
    
%     % only take information corresponding to largest SV
    [U, S, V] = svd(Wd);
    singvals = diag(S);
    fprintf('Removing singular values from %g and smaller\n', singvals(2))
    Wd = U * diag([singvals(1) zeros(1, 2*s_ac.N_b - 1)]) * V';
    
    % check constraint satisfaction 
    assign(W, Wd);
    min(check(Cons_d))
    
    % extract voltage vector and check again
    X = sqrt(diag(Wd)) .* sign(Wd(:,1));
    Wdc =  X * X';
    assign(W,Wdc)
    min(check(Cons_d))
    
    Vdp = X(1:s_ac.N_b) + 1i * X(s_ac.N_b+1:end);
    Vd(union(buses{ii, :}), ii) = Vdp;
end

%% extract voltage vectors for decomposed solutions

Vd = nan(ac.N_b, m);
for ii = 1:m
    
    % get full local matrix 
    Wd = agents(ii).W;
    s_ac = agents(ii).ac;
    W = sdpvar(2*s_ac.N_b);
    Cons_d = [feasibleW(W, wind.P_wf, s_ac), W >= 0];
    
    % only take information corresponding to largest SV
    [U, S, V] = svd(Wd);
    singvals = diag(S);
    fprintf('Removing singular values from %g and smaller\n', singvals(2))
    Wd = U * diag([singvals(1) zeros(1, 2*s_ac.N_b - 1)]) * V';
    
    % check constraint satisfaction 
    assign(W, Wd);
    min(check(Cons_d))
    
    % extract voltage vector and check again
    X = sqrt(diag(Wd)) .* sign(Wd(:,1));
    Wdc =  X * X';
    assign(W,Wdc)
    min(check(Cons_d))
    
    Vdp = X(1:s_ac.N_b) + 1i * X(s_ac.N_b+1:end);
    Vd(union(buses{ii, :}), ii) = Vdp;
end


Vda = nanmean(Vd, 2);
X = [real(Vda); imag(Vda)];
Wd = X * X';

assign(Wf, Wd);
res = check(Cons);

min(res)

%% orrr, take only the own bus set
Vd = nan(ac.N_b, 1);
for ii = 1:m
    
    Wd = agents(ii).broadcast(buses{ii,ii});
    X = sqrt(diag(Wd)) .* sign(Wd(:,1));
    nb = length(X)/2;
    Vd(buses{ii,ii}) = X(1:nb) + 1i*X(nb+1:end);
end

X = [real(Vd); imag(Vd)];
Wd = X * X';

assign(Wf, Wd);
res = check(Cons);

min(res)
