clc
close all
% clear all
yalmip('clear')


%%



ac = load(   'cases30.mat', 'BBUS', 'PGmax', 'QGmax', ...
             'PGmin', 'QGmin', 'C0', 'C1', 'C2', 'PD', 'QD','CG','bl',...
             'bm','Slmax','Vbmax' ,'Vbmin','bg','y','y_','Vbase','Sbase'   ) ;
         

%%
% power network data

NB = 30;   %k
NG = 6;   %g
NL = 41;   %l,m

refIm = 31;
refRe = 1;

% costs

C0 = ac.CG * ac.C0;   % startup cost
C1 = ac.CG * ac.C1 ;  % C1 * PG                                  
C2 = ac.CG * ac.C2 ;    % C2 * PG^2

% power generators limitation

PGmax = ac.CG * ac.PGmax;  %(NB * 1)
QGmax = ac.CG * ac.QGmax;  %(NB * 1)
PGmin = ac.CG * ac.PGmin;  %(NB * 1)
QGmin = ac.CG * ac.QGmin;  %(NB * 1)


% apparent power of demand
 
PD = ac.PD;  %(NB * 1)
QD = ac.QD;  %(NB * 1)


Pinjmax = PGmax - PD;
Pinjmin = PGmin - PD;

Qinjmax = QGmax - QD;
Qinjmin = QGmin - QD;

% compex voltage limitation at every bus

Vbmax = ac.Vbmax; %(NB * 1)
Vbmin = ac.Vbmin;  %(NB * 1)

% lines limitation

Slmax = ac.Slmax;  %(NL * 1)

% admittance of network

Vbase = ac.Vbase;
Sbase = ac.Sbase;
BBUS = ac.BBUS; 

y = ac.y;
y_ = ac.y_ ;

bl = ac.bl;
bm = ac.bm;
bg = ac.bg;


%%
% coding of the buses constraints

for b = 1:NB

    e = zeros(NB,1);
    e(b,1) = 1;  %ek
    
    Y = ( e * e' ) * BBUS;
    
    Mk{b,1} = [ ( e * e' )  , zeros(NB,NB) 
         zeros(NB,NB) ,  ( e * e' )  ];
    
    Yk{b,1} = (1/2) * [ real( Y + Y.' ) , imag( Y.' - Y ) 
                        imag( Y - Y.' ) , real( Y + Y.' ) ];
    
    Yk_{b,1} = (-1/2) * [ imag( Y + Y.' ) , real( Y - Y.' ) 
                          real( Y.' - Y ) , imag( Y + Y.' ) ];
    
   
    if any( b == bg )
        
        A(b,1) = C0(b,1) + ( C1(b,1) * PD(b,1) );
    
        B(b,1) = sqrt( C2(b,1) ) * PD(b,1);
        
    end
    
end

%%
        

for i = 1:NL
    
    l = bl(i);
    m = bm(i);
    
        
    e = zeros(NB,1); %el
    e(l,1) = 1;
    
    ee = zeros(NB,1); %em
    ee(m,1) = 1;
    
    Y = ( ( y_(i,1) + y(i,1) ) * ( e * e' ) ) - ( y(i,1)  * ( e * ee' ) );
        
    Mlm{i,1} = [( e - ee ) * ( e - ee )' , zeros(NB,NB) 
                            zeros(NB,NB) , ( e - ee ) * ( e - ee )'];
            
   Ylm{i,1} = (1/2) * [ real( Y + Y.' ) , imag( Y.' - Y ) 
                        imag( Y - Y.' ) , real( Y + Y.' ) ];
             
    Ylm_{i,1} = (-1/2) * [ imag( Y + Y.' ) , real( Y - Y.' ) 
                           real( Y.' - Y ) , imag( Y + Y.' ) ];
    
    DVlm(i,1) = Vbmax(l,1) + Vbmax(m,1);
 end


 %%
 
 save PreOPF30.mat
 
 
