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

C = [C, W >= 0];
C = [C, W(refIm,refIm) == 0];

%% constraints

for b = 1:NB
       
     C = [C, ( PGmin(b,1) - PD(b,1) ) <= trace( Yk{b,1} * W ) <= ( PGmax(b,1) - PD(b,1) ) ];
                    
     C = [C, ( QGmin(b,1) - QD(b,1) ) <= trace( Yk_{b,1} * W ) <= ( QGmax(b,1) - QD(b,1) ) ];    
     
     C = [C,  ( Vbmin(b,1)^2 ) <= trace( Mk{b,1} * W ) <= ( Vbmax(b,1)^2 ) ];
         
end

%%
 
% for i = 1:NL
%     
%     C = [C, [   ( -( Slmax(i,1) )^2 )        , trace( Ylm{i,1} * W ) ,     trace( Ylm_{i,1} * W ) 
%                 
%                           trace( Ylm{i,1} * W )   ,          -1           ,               0 
%                            
%                           trace( Ylm_{i,1} * W )  ,           0           ,               -1            ] <= 0 ];
%                       
%     C = [C, trace( Mlm{i,1} * W ) <= ( DVlm(i,1)^2 ) ];
%         
% end

%%
% 
% for i = 1:NG
%     
%    b = bg(i,1);
%    
%    C = [C, [ ( C1(b,1) * trace( Yk{b,1} * W ) ) - alfa(b,1) + A(b,1) , ( sqrt( C2(b,1) ) * trace( Yk{b,1} * W ) ) + B(b,1)
%                               
%                   ( sqrt( C2(b,1) ) * trace( Yk{b,1} * W ) ) + B(b,1)     ,                    -1                               ] <= 0 ];
%         
% end
%            
%% Objective

% O = sum(alfa);
O = [];
%%

opt = sdpsettings('solver','sedumi','sedumi.eps',2e-8,'sedumi.bigeps',1e-3,'verbose',1,'shift',0,'relax',0,'savesolverinput',true);

solvesdp(C,O,opt)

%%

Wopt = double(W);
[S,D] = eig(Wopt);
DD = diag(D);
%%%

obj = double(alfa);
OO = sum(obj);
%%%

X=sqrt(diag(Wopt)).*sign(Wopt(:,1));
Xr = X(1:30);
Xi = X(31:60);
%%%

Vbc = Xr + Xi * 1i;
LenVb = abs(Vbc);
AngVb = angle(Vbc) * (180/pi);
%%%

for i = 1:NL
    
    Plm(i,1) = trace( Ylm{i,1} * Wopt );
    Qlm(i,1) = trace( Ylm_{i,1} * Wopt );
    
end
%%%

for b = 1:NB
   
    Pinj(b,1) = trace(Yk{b,1} * Wopt);
    Qinj(b,1) = trace(Yk_{b,1} * Wopt);
        
end
%%%

for i = 1:NG
         
    b = bg(i);
    Pgen(i,1) = ( Pinj(b,1) + PD(b,1) );
    Qgen(i,1) = ( Qinj(b,1) + QD(b,1) );
    Cost(i,1) = ( C2(b,1) * ( Pinj(b,1) + PD(b,1) )^2 ) + ( C1(b,1) * ( Pinj(b,1) + PD(b,1) ) ) + C0(b,1);
    Obj(i,1) = obj(b,1);
    
end
%%%

results30.Wopt = Wopt;
results30.Obj = Obj;
results30.Cost = Cost;
%%%

results30.Xr = Xr;
results30.Xi = Xi;
%%%

results30.LenVb = LenVb;
results30.AngVb = AngVb;
%%%

results30.Pinj = Pinj;
results30.Qinj = Qinj;
%%%

results30.Plm = Plm;
results30.Qlm = Qlm;
%%%

results30.Pgen = Pgen;
results30.Qgen = Qgen;
%%%

save results30.mat results30 
 
 
 
 
