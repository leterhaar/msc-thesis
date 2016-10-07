clc
close all
clear all

%%

mpc = loadcase('case_ieee30');

mpc.branch(:,6) = mpc.branch(:,6) /100 ;

[results, success] = runopf(mpc);

for i=1:30
    
    Xbc(i,1) = results.bus(i,8)*exp(1j*(results.bus(i,9)*(pi/180)));
    Xr(i,1) = real(Xbc(i,1));
    Xi(i,1) = imag(Xbc(i,1));
    
end

resultsM.Cost = results.f;

resultsM.Xr = Xr;
resultsM.Xi = Xi;

resultsM.LenVb = results.bus(:,8);
resultsM.AngVb = results.bus(:,9);

resultsM.Plm = results.branch(:,14);
resultsM.Qlm = results.branch(:,15);

resultsM.Pgen = results.gen(:,2);
resultsM.Qgen = results.gen(:,3);

save resultsMATPower.mat resultsM
 


