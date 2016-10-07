function [sys] = getSystemData(powersystem)
% GETSYSTEMDATA defines network data
% INPUT: powersystem
% 'esp': compact form of the spanish network from "Soler 2010"
% 'rts': IEEE 96 network from "Kirschen 2007"
% 'threbus': 3 Bus network from "Conejo 2009" 
% 'esp_full': spanish network from "Soler 2010"
% OUTPUT: struct with network data

switch lower(powersystem)
    case 'esp'
        %% 'esp': compact form of the spanish network from "Soler 2010"
        % Indices and numbers
        sys.N_B = 6; % buses n=1...N_B
        sys.N_G = 7; % generators i=1...N_G
        sys.N_L = 4; % loads j=1...N_L
        sys.N_l = 10; % lines nl=1...N_l

        % M_WP: 1 x bus, 1's where WP infeed happens
        sys.M_WP =zeros(sys.N_B,1);
        sys.M_WP(3) = 1;
        
        % Production limits and loads
        sys.Pi_min = diag([75 75 75 5 10 3 15]);
        sys.Pi_max = diag([800 600 200 100 180 160 180]);
                
        sys.load_prof_day = [1194,1162,1141,1135,1152,1268,1433,1542,1595,1642,1665,1657,...
                             1606,1535,1537,1555,1575,1708,1750,1717,1594,1568,1385,1296];
        
        sys.load_perc = [0.14 0.51 0.26 0.09];

        % Price values
        sys.lambda_itSU = diag([3606 4207 4508 240 361 421 0]);
        sys.lambda_itG = [17 16 18 20 19 19 16]';
        sys.lambda_jtL = 0*ones(sys.N_L,1);
        sys.V_jtLOL = 1000;
        sys.V_spill = 0;
        sys.V_spill_slack = 1000;
        
        sys.C_itRU = round([1/6*sys.lambda_itG(1:3);2/5*sys.lambda_itG(4:7)]*10)/10;
        sys.C_itRD = round([1/6*sys.lambda_itG(1:3);2/5*sys.lambda_itG(4:7)]*10)/10;
        sys.C_itRNS = round([1/7*sys.lambda_itG(1:3);1/3*sys.lambda_itG(4:7)]*10)/10;
        sys.C_jtRU = [70*ones(1,sys.N_L)]';
        sys.C_jtRD = [70*ones(1,sys.N_L)]';
        
        % Network data:
        % M_G: bus x generators, 1's where a generator is connected
        sys.M_G = zeros(sys.N_B,sys.N_G);
        sys.M_G(1,1) = 1; % bus 1
        sys.M_G(2,2:4) = 1; % bus 2
        sys.M_G(3,5) = 1; % bus 3
        sys.M_G(4,6) = 1; % bus 4
        sys.M_G(6,7) = 1; % bus 6

        % M_L: bus x loads, 1's where a load is connected
        sys.M_L = zeros(sys.N_B,sys.N_L);
        sys.M_L(2,1) = 1; % bus 2
        sys.M_L(3,2) = 1; % bus 3
        sys.M_L(4,3) = 1; % bus 4
        sys.M_L(5,4) = 1; % bus 5
        
        % The connectivity matrix (lines x buses)
        sys.M_c = zeros(sys.N_l,sys.N_B);
        sys.M_c(1:3,1) = 1; % from bus 1
        sys.M_c([4:5,9:10],2) = 1; % from bus 2
        sys.M_c(6:7,3) = 1; % from bus 3
        sys.M_c(8,4) = 1; % from bus 4
        sys.M_c(1,2) = -1; % to bus 2
        sys.M_c(2:3,3) = -1; % to bus 3
        sys.M_c(4:6,4) = -1; % to bus 4
        sys.M_c(7:8,5) = -1; % to bus 5
        sys.M_c(9:10,6) = -1; % to bus 6

        % V_base [kV]:
        sys.V_base = 400;
        
        % Imaginary part of the line admittances
        X = [8.1 9.5 9.5 10.2 10.2 7.6 4.8 4.8 6 6];
        sys.Y_im = diag((pi/180)*sys.V_base^2./X);


        % Transmission capacity limits
        Al=[266.5 799.5 799.5 533 533 266.5 266.5 266.5 160 160];
%       Bl=[183.5 550.5 550.5 367 367 183.5 183.5 183.5 110 110];
        sys.f_lineMax = Al'; % line 1-10
        
    case 'rts'
        %% 'rts': IEEE 96 network from "Kirschen 2007"
        % Indices and numbers
        sys.N_B = 24; % buses n=1...N_B
        sys.N_G = 26; % generators i=1...N_G
        sys.N_L = 17; % loads j=1...N_L
        sys.N_l = 38; % lines nl=1...N_l

        % M_WP: 1 x bus, 1's where WP infeed happens
        sys.M_WP =zeros(sys.N_B,1);
        sys.M_WP(16) = 1;
        
        % Production limits and loads
        sys.Pi_min = diag([2.4,2.4,2.4,2.4,2.4,4,4,4,4,15.2,15.2,15.2,15.2,2.5,2.5,2.5,54.25,54.25,54.25,54.25,68.95,68.95,68.95,140,100,100]);
        sys.Pi_max = diag([12,12,12,12,12,20,20,20,20,76,76,76,76,100,100,100,155,155,155,155,197,197,197,350,400,400]);
        
        sys.load_prof_day = [1700,1730,1690,1700,1750,1850,2000,2430,2540,2600,2670,2590,...
                         2590,2550,2620,2650,2550,2530,2500,2550,2600,2480,2200,1840];
%         load_prof_day = [1430,1450,1400,1350,1350,1470,1710,2060,2300,2380,2290,2370,...
%                            2290,2260,2190,2130,2190,2200,2300,2340,2300,2180,1910,1650];
        sys.load_perc = (1/100)*[3.8,3.4,6.3,2.6,2.5,4.8,4.4,6.0,6.1,6.8,...% Load 1-10
                                 9.3,6.8,11.1,3.5,11.7,6.4,4.5]; % Load 11-17

        % Price values
        sys.lambda_itSU = diag([68 68 68 68 68 5 5 5 5 655.6 655.6 655.6 655.6 566 566 566 1048.3 1048.3 1048.3 1048.3 775 775 775 4468 5000 5000]);
        sys.lambda_itG = [25.5472,25.6753,25.8027,25.9318,26.0611,37.551,37.6637,37.777,37.8896,13.3272,13.3538,13.3805,13.4073,18,18.1,18.2,10.694,10.7154,10.7367,10.7583,23,23.1,23.2,10.8616,7.4921,7.5031]';
        sys.lambda_jtL = 0*ones(sys.N_L,1);
        sys.V_jtLOL = 1000;
        sys.V_spill = 0;
        sys.V_spill_slack = 1000;
        
        sys.C_itRU = round(10*[2/5*sys.lambda_itG(1:9); 1/6*sys.lambda_itG(10:end)])/10;
        sys.C_itRD = round(10*[2/5*sys.lambda_itG(1:9); 1/6*sys.lambda_itG(10:end)])/10;
        sys.C_itRNS = round(10*[1/3*sys.lambda_itG(1:9); 1/7*sys.lambda_itG(10:end)])/10;
        sys.C_jtRU = [70*ones(1,sys.N_L)]';
        sys.C_jtRD = [70*ones(1,sys.N_L)]';
        
        % Network data:
        % M_G: bus x generators, 1's where a generator is connected
        sys.M_G = zeros(sys.N_B,sys.N_G);
        g_at_bus = [15 15 15 15 15 1 1 2 2 1 1 2 2 7 7 7 15 16 23 23 13 13 13 23 18 21];
        for g = 1:sys.N_G        
            sys.M_G(g_at_bus(g),g) = 1; % connect generator g
        end
        
        % M_L: bus x loads, 1's where a load is connected
        sys.M_L = zeros(sys.N_B,sys.N_L);
        L_at_bus = [1:10 13:16 18:20];
        for L = 1:sys.N_L       
            sys.M_L(L_at_bus(L),L) = 1; % connect load L
        end
        
        % The connectivity matrix (lines x buses)
        sys.M_c = zeros(sys.N_l,sys.N_B);
        from_bus = [1,1,1,2,2,3,3,4,5,6,7,8,8,9,9,10,10,11,11,12,12,13,14,15,15,15,15,16,16,17,17,18,18,19,19,20,20,21];
        to_bus = [2,3,5,4,6,9,24,9,10,10,8,9,10,11,12,11,12,13,14,13,23,23,16,16,21,21,24,17,19,18,22,21,21,20,20,23,23,22];
        for l=1:sys.N_l
            sys.M_c(l,from_bus(l)) = 1; %line from bus in left Networkside
            sys.M_c(l,to_bus(l)) = -1; %line to bus in left Networkside
        end
   
        % Imaginary part of the line admittances
        X_pu = [0.0139000000000000,0.211200000000000,0.0845000000000000,0.126700000000000,0.192000000000000,0.119000000000000,0.0839000000000000,0.103700000000000,0.0883000000000000,0.0605000000000000,0.0614000000000000,0.165100000000000,0.165100000000000,0.0839000000000000,0.0839000000000000,0.0839000000000000,0.0839000000000000,0.0476000000000000,0.0418000000000000,0.0476000000000000,0.0966000000000000,0.0865000000000000,0.0389000000000000,0.0173000000000000,0.0490000000000000,0.0490000000000000,0.0519000000000000,0.0259000000000000,0.0231000000000000,0.0144000000000000,0.105300000000000,0.0259000000000000,0.0259000000000000,0.0396000000000000,0.0396000000000000,0.0216000000000000,0.0216000000000000,0.0678000000000000];
        sys.S_base = 100;% [MVA]
        sys.V_base = 138;% [kV]
%         % V_base [kV] diffrent for each bus:
%         sys.V_base = 230*ones(1,sys.N_B);%[138,138,138,138,138,138,138,138,138,138,230,230,230,230,230,230,230,230,230,230,230,230,230,230];
%         M_cp = sys.M_c;
%         M_cp(M_cp<0) = 0; %+1 where M_c is 1, else 0
%         M_cn = -sys.M_c;
%         M_cn(M_cn<0) = 0; %+1 where M_c is -1, else 0        
%         X = 1/sys.S_base * (M_cp*sys.V_base') .*  (M_cn*sys.V_base') .* X_pu'; 
        X = X_pu*sys.V_base^2/sys.S_base;
        sys.Y_im = diag((pi/180)*sys.V_base^2./X);
        
        % Transmission capacity limits
        Al = [175,175,175,175,175,175,400,175,175,175,175,175,175,400,400,400,400,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500];
        sys.f_lineMax = Al'; % line 1-N_l
        
    case 'threebus'
        %% 'threbus': 3 Bus network from "Conejo 2009" 
        % Indices and numbers
        sys.N_B = 3; % buses n=1...N_B
        sys.N_G = 3; % generators i=1...N_G
        sys.N_L = 2; % loads j=1...N_L
        sys.N_l = 3; % lines nl=1...N_l

        % M_WP: 1 x bus, 1's where WP infeed happens
        sys.M_WP =zeros(sys.N_B,1);
        sys.M_WP(2) = 1;
        
        % Production limits and loads
        sys.Pi_min = diag([10 10 10]);
        sys.Pi_max = diag([220 180 300]);
                
        sys.load_prof_day = 170*ones(1,24);
        
        sys.load_perc = [2 1];

        % Price values
        sys.lambda_itSU =  diag([100 100 100]);
        sys.lambda_itG = [30 40 50]';
        sys.lambda_jtL = 0*ones(sys.N_L,1);
        sys.V_jtLOL = 1000;
        sys.V_spill = 0;
        sys.V_spill_slack = 1000;
        sys.C_itRU = [5 7 8]';
        sys.C_itRD = [5 7 8]';
        sys.C_itRNS = [4.5 5.5 7]';
        sys.C_jtRU = [70 70]';
        sys.C_jtRD = [70 70]';
        
        % Network data, consider the 3-bus system of Fig.1
        % M_G: bus x generators, 1's where a generator is connected
        sys.M_G = eye(sys.N_G);
        
        % M_L: bus x loads, 1's where a load is connected
        sys.M_L = [0 0;1 0;0 1];
        
        % The connectivity matrix (lines x buses)
        sys.M_c = [1 -1 0; 1 0 -1; 0 1 -1];
        
        % V_base [kV] diffrent for each bus:
        sys.Vbase = 120;
        sys.Sbase = 200;
        
        % Imaginary part of the line admittances
        sys.X = 0.13; % Ohm
        sys.Y_im = (sys.Vbase^2/sys.X)*eye(sys.N_l); % (pi/180)*(sys.Vbase^2/sys.X)*eye(sys.N_l);
        % Transmission capacity limits
        sys.f_lineMax = [110; 110; 110];
    
    case 'spanish_full'
        %% 'esp_full': spanish network from "Soler 2010"
        % Indices and numbers
        sys.N_B = 6; % buses n=1...N_B
        sys.N_G = 18; % generators i=1...N_G
        sys.N_L = 4; % loads j=1...N_L
        sys.N_l = 10; % lines nl=1...N_l

        % M_WP: 1 x bus, 1's where WP infeed happens
        sys.M_WP =zeros(sys.N_B,1);
        sys.M_WP(3) = 1;
        
        % Production limits and loads
        sys.Pi_min = diag([75*ones(1,8),5*ones(1,2) 10*ones(1,3) 3*ones(1,4) 15]);
        sys.Pi_max = diag([170 190 210 230 180 200 220 200 50*ones(1,2) 60*ones(1,3) 40*ones(1,4) 180]);
        
        sys.load_prof_day = [1194,1162,1141,1135,1152,1268,1433,1542,1595,1642,1665,1657,...
                             1606,1535,1537,1555,1575,1708,1750,1717,1594,1568,1385,1296];
        
        sys.load_perc = [0.14 0.51 0.26 0.09];

        % Price values
        sys.lambda_itSU = diag([3606*ones(1,4) 4207*ones(1,3) 4508 240*ones(1,2) 361*ones(1,3) 421*ones(1,4) 0]);
        sys.lambda_itG = [17*ones(1,4) 16*ones(1,3) 18 20*ones(1,2) 19*ones(1,3) 19*ones(1,4) 16]';
        sys.lambda_jtL = 0*ones(sys.N_L,1);
        sys.V_jtLOL = 1000;
        sys.V_spill = 0;
        sys.V_spill_slack = 1000;
        
        sys.C_itRU = round([1/6*sys.lambda_itG(1:8);2/5*sys.lambda_itG(9:18)]*10)/10;
        sys.C_itRD = round([1/6*sys.lambda_itG(1:8);2/5*sys.lambda_itG(9:18)]*10)/10;
        sys.C_itRNS = round([1/7*sys.lambda_itG(1:8);1/3*sys.lambda_itG(9:18)]*10)/10;
        sys.C_jtRU = [70*ones(1,sys.N_L)]';
        sys.C_jtRD = [70*ones(1,sys.N_L)]';

        % Network data:
        % M_G: bus x generators, 1's where a generator is connected
        sys.M_G = zeros(sys.N_B,sys.N_G);
        sys.M_G(1,1:4) = 1; % bus 1
        sys.M_G(2,5:10) = 1; % bus 2
        sys.M_G(3,11:13) = 1; % bus 3
        sys.M_G(4,14:17) = 1; % bus 4
        sys.M_G(6,18) = 1; % bus 6

        % M_L: bus x loads, 1's where a load is connected
        sys.M_L = zeros(sys.N_B,sys.N_L);
        sys.M_L(2,1) = 1; % bus 2
        sys.M_L(3,2) = 1; % bus 3
        sys.M_L(4,3) = 1; % bus 4
        sys.M_L(5,4) = 1; % bus 5
        
        % The connectivity matrix (lines x buses)
        sys.M_c = zeros(sys.N_l,sys.N_B);
        sys.M_c(1:3,1) = 1; % from bus 1
        sys.M_c([4:5,9:10],2) = 1; % from bus 2
        sys.M_c(6:7,3) = 1; % from bus 3
        sys.M_c(8,4) = 1; % from bus 4
        sys.M_c(1,2) = -1; % to bus 2
        sys.M_c(2:3,3) = -1; % to bus 3
        sys.M_c(4:6,4) = -1; % to bus 4
        sys.M_c(7:8,5) = -1; % to bus 5
        sys.M_c(9:10,6) = -1; % to bus 6

        % V_base [kV]
        sys.V_base = 400;
        
        % Imaginary part of the line admittances
        X = [8.1 9.5 9.5 10.2 10.2 7.6 4.8 4.8 6 6];
        sys.Y_im = diag((pi/180)*sys.V_base^2./X);
        % Transmission capacity limits
        Al=[266.5 799.5 799.5 533 533 266.5 266.5 266.5 160 160];
%       Bl=[183.5 550.5 550.5 367 367 183.5 183.5 183.5 110 110];
        sys.f_lineMax = Al'; % line 1-10
    otherwise
            disp('Not a valid case! in getSystemData.m');
end
