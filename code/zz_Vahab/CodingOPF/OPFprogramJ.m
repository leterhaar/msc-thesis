clc
close all
clear all
yalmip('clear')

%%

load PreOPF30.mat

%%
% define variables

alf = sdpvar(NG,1,'full');
alfa = ac.CG * alf;

W = sdpvar(2*NB,2*NB,'symmetric');

%%

C = [];

C = C + set ( W >= 0 );

%% constraints

for b = 1:NB
    
    
     C = C + set( ( PGmin(b,1) - PD(b,1) ) <= trace( Yk{b,1} * W ) <= ( PGmax(b,1) - PD(b,1) ) );
                    
     C = C + set( ( QGmin(b,1) - QD(b,1) ) <= trace( Yk_{b,1} * W ) <= ( QGmax(b,1) - QD(b,1) ) );    
     
     C = C + set(  ( Vbmin(b,1)^2 ) <= trace( Mk{b,1} * W ) <= ( Vbmax(b,1)^2 ) );
         
end

%%
        

for i = 1:NL
    
    C = C + set( [   ( -( Slmax(i,1) )^2 )        , trace( Ylm{i,1} * W ) ,     trace( Ylm_{i,1} * W ) 
                
                          trace( Ylm{i,1} * W )   ,          -1           ,               0 
                           
                          trace( Ylm_{i,1} * W )  ,           0           ,               -1            ] <= 0 );

    C = C + set( trace( Mlm{i,1} * W ) <= ( DVlm(i,1)^2 ) );
        
end

%%

for i = 1:NG
    
   b = bg(i,1);
   
   C = C + set( [   ( C1(b,1) * trace( Yk{b,1} * W ) ) - alfa(b,1) + A(b,1) , ( sqrt( C2(b,1) ) * trace( Yk{b,1} * W ) ) + B(b,1)
                              
                    ( sqrt( C2(b,1) ) * trace( Yk{b,1} * W ) ) + B(b,1)     ,                    -1                               ] <= 0 );
         
   
        
end
           
%% Objective

O = sum(alfa);

%%

opt = sdpsettings('solver','sedumi','sedumi.eps',5e-6,'sedumi.bigeps',1e-3,'verbose',1,'shift',0,'relax',0);

solvesdp(C,O,opt)

%%

obj = double(alfa);
OO = sum(obj);
%%%

W1 = double(W);
[S,D] = eig(W1);
DD = diag(D);

p1 = DD(59,1);
p2 = DD(60,1);
E = S(:,59);

Wopt = (p1+p2)*(E*E');
%%%

X = sqrt(diag(Wopt)).*sign(Wopt(:,1));
Xr = X(1:30);
Xi = X(31:60);
%%%

Vbc = Xr + Xi * 1i;

LenVb = abs(Vbc);
AngVb = angle(Vbc) * (180/pi);
%%%

for i = 1:41
    
    Plm(i,1) = trace( Ylm{i,1} * Wopt );
    Qlm(i,1) = trace( Ylm_{i,1} * Wopt );
    
end
%%%

for b = 1:30
   
    Pinj(b,1) = trace(Yk{b,1} * W1);
    Qinj(b,1) = trace(Yk_{b,1} * W1);
        
end
%%%

for i = 1:6
         
    b = bg(i);
    Pgen(i,1) = ( Pinj(b,1) + PD(b,1) );
    Qgen(i,1) = ( Qinj(b,1) + QD(b,1) );
    Cost(i,1) = ( C2(b,1) * ( Pinj(b,1) + PD(b,1) )^2 ) + ( C1(b,1) * ( Pinj(b,1) + PD(b,1) ) ) + C0(b,1);
    Obj(i,1) = obj(b,1);
    
end
%%%

resultsJ.Wopt = Wopt;
%%%

resultsJ.Obj = Obj;
resultsJ.Cost = Cost;
%%%

resultsJ.Xr = Xr;
resultsJ.Xi = Xi;
%%%

resultsJ.LenVb = LenVb;
resultsJ.AngVb = AngVb;
%%%

resultsJ.Pinj = Pinj;
resultsJ.Qinj = Qinj;
%%%

resultsJ.Plm = Plm;
resultsJ.Qlm = Qlm;
%%%

resultsJ.Pgen = Pgen;
resultsJ.Qgen = Qgen;
%%%

save resultsJ.mat resultsJ 
