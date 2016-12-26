%% test SDPs for ACCA algorithm

% start freshh
clear 
yalmip('clear');
clc

% add paths 
if not(exist('read_sdpa.m', 'file'))
    addpath('../misc');
    addpath('../sdplib');
end

sdp = read_sdpa('hinf1.dat-s');

obj = @(X) trace(sdp.C * X);
C = [];
for i = 1:length(sdp.As)
    C = [C, trace(sdp.As{i} * sdp.X) == sdp.bs(i)];
end
status = optimize([C, sdp.X >= 0], obj(sdp.X), ...
                                    sdpsettings('solver', 'sedumi'));
verify(not(status.problem), status.info);