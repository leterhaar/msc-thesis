clc
clear all
close all

%%

load('results30');
load('resultsJ');
load('resultsMATpower');

set(0,'defaultaxesfontsize',16);

mycolors = [1 0 0; 0 1 0; 0 0 1];

%%

TotalCost = [resultsM.Cost  sum(results30.Cost)  sum(resultsJ.Cost)]/1000;

Pgen = [resultsM.Pgen/100 results30.Pgen resultsJ.Pgen]*100;
Qgen = [resultsM.Qgen/100 results30.Qgen resultsJ.Qgen]*100;

Plm = [resultsM.Plm/100 results30.Plm resultsJ.Plm]*100;
Qlm = [resultsM.Qlm/100 results30.Qlm resultsJ.Qlm]*100;

LenVb = [resultsM.LenVb results30.LenVb resultsJ.LenVb];
AngVb = [resultsM.AngVb results30.AngVb resultsJ.AngVb];

%%

figure(1)
bar(TotalCost,'LineWidth',2,'LineStyle',':')
set(gca,'XTickLabel',{'MatPower';'Opt.4';'Opt.3'})
% colormap summer % Change the color scheme
colormap(mycolors);
grid on 
ylabel('[$/MW]')
ylim([8.5 9.5])
title('Total Cost of OPF problem');

%%

figure(2)
subplot(2,1,1)
bar(Pgen,'LineWidth',2,'LineStyle',':')
set(gca,'XTickLabel',{'Gen1';'Gen2';'Gen3';'Gen4';'Gen5';'Gen6'})
%colormap summer % Change the color scheme
colormap(mycolors);
% legend('MatPower','Opt.4','Opt.3')
grid on 
ylabel('[MW]')
title('Active Power of Generate Units');

%%

figure(2)
subplot(2,1,2)
bar(Qgen,'LineWidth',2,'LineStyle',':')
set(gca,'XTickLabel',{'Gen1';'Gen2';'Gen3';'Gen4';'Gen5';'Gen6'})
% colormap summer % Change the color scheme
colormap(mycolors);
% legend('MatPower','ImprovedSDP OPF','SDP OPF')
grid on 
ylabel('[MW]')
ylim([0 50])
title('Reactive Power of Generate Units');

%%

figure(3)
subplot(2,1,1)
bar(Plm)
% colormap summer % Change the color scheme
colormap(mycolors);
legend('MatPower','Opt.4','Opt.3')
grid on 
ylabel('[MW]')
xlabel('Number of Lines')
xlim([0 42])
ylim([0 145])
title('Active Power flow in each line');

%%

figure(3)
subplot(2,1,2)
bar(Qlm)
% colormap summer % Change the color scheme
colormap(mycolors);
% legend('MatPower','ImprovedSDP OPF','SDP OPF')
grid on 
ylabel('[MW]')
xlabel('Number of Lines')
xlim([0 42])
title('Reactive Power flow in each line');

%%

figure(4)
subplot(2,1,1)
bar(LenVb)
% colormap summer % Change the color scheme
colormap(mycolors);
% legend('MatPower','ImprovedSDP OPF','SDP OPF')
grid on 
ylabel('Kvolt')
xlabel('Number of Buses')
xlim([0 31])
ylim([0.8 1.1])
title('Voltage Length in every Buses');

%%

figure(4)
subplot(2,1,2)
bar(AngVb)
% colormap summer % Change the color scheme
colormap(mycolors);
legend('MatPower','Opt.4','Opt.3')
grid on 
ylabel('Degree')
xlabel('Number of Buses')
xlim([0 31])
% ylim([125 150])
title('Voltage Angle in every Buses');

