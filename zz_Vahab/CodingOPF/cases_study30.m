clc
clear all
close all

%%

mpc = loadcase('case_ieee30');

Vbase = mpc.bus(:,10);
Sbase = mpc.baseMVA;
Vbase = Vbase(1,1);

%% base data

Vbmax = mpc.bus(:,12) ;  %[p.u]
Vbmin = mpc.bus(:,13) ;  %[p.u]

PD = mpc.bus(:,3) / Sbase;   %[MW]
QD = mpc.bus(:,4)/ Sbase;   %[MVAr]

gsh = mpc.bus(:,5)/ Sbase;
bsh = mpc.bus(:,6)/ Sbase;
ysh = gsh + bsh *1i;

%%

bg = mpc.gen(:,1);

CG = zeros(30,6);

for c = 1:30
    for cc = 1:length(bg)
        if c == bg(cc,1)
            CG(c,cc) = 1;
        end
    end
end


PGmax = mpc.gen(:,9) / Sbase;   %[MW]
PGmin = mpc.gen(:,10) / Sbase;   %[MW]


QGmax = mpc.gen(:,4) / Sbase;   %[MVAr]
QGmin = mpc.gen(:,5) / Sbase;   %[MVAr]


C0 = mpc.gencost(:,7) * Sbase;   %[$/H]
C1 = mpc.gencost(:,6) * Sbase;   %[$/H]
C2 = mpc.gencost(:,5) * Sbase^2;   %[$/H]


%%

Slmax = mpc.branch(:,6) /Sbase^2 ;

bl = mpc.branch(:,1);
bm = mpc.branch(:,2);



x = ( mpc.branch(:,3) ) + mpc.branch(:,4) * 1i;


for k = 1:length(x)
    
    y(k,1) = ( 1/x(k,1) ) ;   % admittance 
    yy(k,1) = (1/2) * mpc.branch(k,5)  * 1i;   %shunt
end
    
for s =1:length(x)
    
    sh = bl(s);
     y_(s,1) = yy(sh,1) + ysh(sh,1);
     
end

%%

[YBUS, YF, YT] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);

BF = full(YF);
BT = full(YT);
BBUS = full(YBUS) ;

%%

[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);

Bf = full(Bf);
Bbus = full(Bbus) ;

%%

caseformat





save cases30.mat Bf Bbus PGmax PGmin QGmax QGmin C0 C1 C2 PD QD Vbmax Vbmin BF BT BBUS CG bl bm Slmax bg Vbase Sbase y y_













