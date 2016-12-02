clear
clc
yalmip('clear');
reps = 50;
[times(1:reps).direct] = deal(nan);
[times(1:reps).function] = deal(nan);

d = 5;
lp.x = sdpvar(d,1);
lp.z = randn(d,1);
a = 0.2;
lp.f = randn(1,d) * lp.x;
opt_settings = sdpsettings('verbose', 0);
p = progress('Timing LP application', reps);
for i = 1:reps
    % time direct application
    tic
    f_prox = lp.f + 1/(2*a)*norm(lp.x - lp.z)^2;
    info = optimize([], f_prox, opt_settings);
    assert(not(info.problem), 'Problem optimizing: %s', info.info);
    x_direct = value(lp.x);
    times(i).direct = toc;
    clear f_prox
    
    % time proximal application
    tic
    x_function = prox(lp.f, a, lp.z, lp.x);
    times(i).function = toc;
    
    assert(all_close(x_direct, x_function), 'Not close');
    
    p.ping();
end
%%
avg_direct = mean([times.direct]);
avg_function = mean([times.function]);
fprintf(['Average execution time (%i reps): \n'...
         'Direct \t\t%5g sec (%2g%%) \nFunctional \t%5g sec (%.3g%%)\n'], ...
         reps, avg_direct, 100, avg_function,  avg_function/avg_direct * 100);

%%
yalmip('clear');
[sdp_times(1:reps).direct] = deal(nan);
[sdp_times(1:reps).function] = deal(nan);
sdp = read_sdpa('../sdplib/hinf1.dat-s');

z = rand(size(sdp.X,1),1);
sdp.Z = z*z';
sdp.f = trace(sdp.C*sdp.X);

p = progress('Timing SDP application', reps);
for i = 1:reps
    
    % time direct application
    tic
    f_prox = sdp.f + 1/(2*a)*norm(sdp.X - sdp.Z, 'fro')^2;
    info = optimize([], f_prox, opt_settings);
    assert(not(info.problem), 'Problem optimizing: %s', info.info);
    x_direct = value(sdp.X);
    sdp_times(i).direct = toc;
    clear f_prox
    
    % time proximal application
    tic
    x_function = prox(sdp.f, a, sdp.Z, sdp.X);
    sdp_times(i).function = toc;
    
    assert(all_close(x_direct, x_function), 'Not close');
    
    p.ping();
end
%%
avg_direct = mean([sdp_times.direct]);
avg_function = mean([sdp_times.function]);
fprintf(['Average execution time (%i reps): \n'...
         'Direct \t\t%5g sec (%2g%%) \nFunctional \t%5g sec (%.3g%%)\n'], ...
         reps, avg_direct, 100, avg_function,  avg_function/avg_direct * 100);
     
%% plot
initfig('Timing');
plot([times.direct], '-s', 'LineWidth', 2, 'color', orange);
plot([times.function], '-s', 'LineWidth', 2, 'color', blue);
plot([sdp_times.direct], '-.o', 'LineWidth', 2, 'color', orange);
plot([sdp_times.function], '-.o', 'LineWidth', 2, 'color', blue);
    
    