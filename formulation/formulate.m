classdef formulate < handle
    %% Formulate the problem 
    % 
    % This class handles the construction of constraints and objective
    % functions in preparation for the optimization. After the solution has
    % been found, the result can be checked using the exact same functions
    % to validate the constraint satisfaction.
    
    properties 
        name;               % formulation name
        checked_constraints;% Struct with all the constraints (from evaluate)
        details;            % string with details of evaluation
    end
    
    properties (Transient = true)
        prepped_constraints;% Struct with all the constraints (from prepare, SDP vars)
        Constraints;        % LMIs with all the constraints (SDP vars)
        Objective;          % Objective function (SDP vars)
    end

    methods
        
        function obj = formulate(name)
            %% Define the type of formulation
            %
            % Arguments 
            % ---------
            % name : string, name of the formulation. Default 'AC_OPF_RS_affine'
            if nargin < 1
                name = 'P3';
            end
            
            obj.name = name;
        end
        
        function cons = get_constraints(obj, varargin)
            %% Return a struct containing the fixed constraints
            %
            % Arguments
            % ----------
            % see help of prepare() or evaluate()
            %
            % Returns
            % -------
            % cons : a struct with the following cell fields:
            %   bus : bus number corresponding to that constraint
            %   scenario : scenario             "           ", 0=forecast
            %   description : string with a description of the constraint
            %   type : 1) box or 2) psd or 3) equality
            %   lower : lower limit, empty for psd constraint
            %   middle : evaluation of middle part
            %   upper : upper limit, empty for psd constraint
            addpath('../misc');

            if strcmpi(obj.name, 'P1')
            %% P1
                varargin = varargin{1};
                ac = varargin{1};   
                wind = varargin{2};
                t = varargin{3};
                decision_vars = varargin{4};                
                
                % decision variables
                W_f = decision_vars{1};
                W_s = decision_vars{2};
                R_us = decision_vars{3};
                R_ds = decision_vars{4};
                d_us = decision_vars{5};
                d_ds = decision_vars{6};

                % determine number of constraints and preallocate
                [~, N] = size(wind.P_w);
                N_c = 3*ac.N_b + N*(ac.N_b*3 + 2*ac.N_G + 2) + 4;
                [cons(1:N_c).bus] = deal(0); % preallocate struct
                [cons(1:N_c).scenario] = deal(0);
                k_ref = ac.N_b + ac.refbus;
                
                j = 0;
                for k = 1:ac.N_b

                    % (1) P_inj
                    j = j + 1;
                    cons(j).lower =     ac.P_min(k);  
                    cons(j).middle =    trace(ac.Y_k(k)*W_f)       ...
                                        + ac.P_D(t, k)                ...
                                        - ac.C_w(k)*wind.P_wf(t);
                    cons(j).upper =     ac.P_max(k);
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(1) P_inj';
                    
                    % (2) Q_inj
                    j = j+1;
                    cons(j).lower =     ac.Q_min(k);  
                    cons(j).middle =    trace(ac.Ybar_k(k)*W_f)    ...
                                        + ac.Q_D(t, k);
                    cons(j).upper =     ac.Q_max(k);
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(2) Q_inj';
                    
                    % (3) V_bus
                    j = j+1;
                    cons(j).lower =     ac.V_min(k)^2;  
                    cons(j).middle =    trace(ac.M_k(k)*W_f);
                    cons(j).upper =     ac.V_max(k)^2;
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(3) V_bus';
                    
                    
                    
                end
                
                for i = 1:N
                
                    for k = 1:ac.N_b

                        % (1") P_inj scen
                        j = j + 1;
                        cons(j).lower =     ac.P_min(k);  
                        cons(j).middle =    trace(ac.Y_k(k)*W_s(:,:,i)) ...
                                            + ac.P_D(t, k)                 ...
                                            - ac.C_w(k)*wind.P_w(t, i);
                        cons(j).upper =     ac.P_max(k);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(1") P_inj';
                                           
                        % (2") Q_inj scen
                        j = j+1;
                        cons(j).lower =     ac.Q_min(k);  
                        cons(j).middle =    trace(ac.Ybar_k(k)*W_s(:,:,i)) ...
                                            + ac.Q_D(t, k);
                        cons(j).upper =     ac.Q_max(k);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(2") Q_inj';
                        
                        % (3") V_bus scen
                        j = j+1;
                        cons(j).lower =     ac.V_min(k)^2;  
                        cons(j).middle =    trace(ac.M_k(k)*W_s(:,:,i));
                        cons(j).upper =     ac.V_max(k)^2;
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(3") V_bus';
                        
                        
                                              
                    end
                    
                    for l = 1:ac.N_G
                        
                        % reserve power is always between R_us and R_ds
                        j = j + 1;
                        k = ac.Gens(l);
                        cons(j).lower =     -1*R_ds(l);
                        cons(j).middle =    trace(ac.Y_k(k)*(W_s(:,:,i)-W_f)) ...
                                             - ac.C_w(k)*wind.P_m(t, i);
                        cons(j).upper =     R_us(l);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(7) bound R';
                        
                        % relate W_s and W_f through d_ds and d_us
                        j = j + 1;
                        cons(j).middle =    trace(ac.Y_k(k)*(W_s(:,:,i)-W_f)) ...
                                            - ac.C_w(k)*wind.P_m(t, i) - ...
                                            (d_us(l) * max(0, -wind.P_m(t, i)) ...
                                            + d_ds(l) * max(0, wind.P_m(t, i)));
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 3;
                        cons(j).description = '(8) Relation W_s and d';
 
                    end
                    
                    
                    % refbus angle on W_s
                    j = j + 1;                  
                    cons(j).middle = W_s(k_ref, k_ref, i);
                    cons(j).description = '(5") Refbus angle W_s';
                    cons(j).type = 3;
                    cons(j).scenario = i;
                    
                    % PSD on W_s
                    j = j + 1;
                    cons(j).middle = W_s(:,:,i);
                    cons(j).description = '(6") PSD-ness W_s';
                    cons(j).type = 2;
                    cons(j).scenario = i;
                    
                end
                
                % PSD constraints
                j = j + 1;
                cons(j).middle = W_f;
                cons(j).description = '(6) PSD-ness W_f';
                cons(j).type = 2;
                cons(j).scenario = 0;
                               
                % reference bus angle constraints
                j = j + 1;
                cons(j).middle = W_f(k_ref, k_ref);
                cons(j).description = '(5) Refbus angle W_f';
                cons(j).type = 3;
                cons(j).scenario = 0;
                
                % nonnegativity constraints on R limits
                j = j + 1;
                cons(j).middle = R_us;
                cons(j).description = '(9a) NN R_up';
                cons(j).type = 2;
                cons(j).scenario = 0;
                j = j + 1;
                cons(j).middle = R_ds;
                cons(j).description = '(9b) NN R_down';
                cons(j).type = 2;
                cons(j).scenario = 0;
                
                
                
                
            elseif strcmpi(obj.name, 'P3')
            %% P3
                
                varargin = varargin{1};
                ac = varargin{1};   
                wind = varargin{2};
                t = varargin{3};
                decision_vars = varargin{4};

                % decision variables
                W_f = decision_vars{1};
                W_m = decision_vars{2};
                R_us = decision_vars{3};
                R_ds = decision_vars{4};


                % determine number of constraints and preallocate
                [~, N] = size(wind.P_w);
                N_c = 3*ac.N_b + N*(ac.N_b*3 + ac.N_G + 2) + 4;
                [cons(1:N_c).bus] = deal(0); % preallocate struct
                [cons(1:N_c).scenario] = deal(0);
                k_ref = ac.N_b + ac.refbus;
                
                j = 0;
                for k = 1:ac.N_b

                    % (1) P_inj
                    j = j + 1;
                    cons(j).lower =     ac.P_min(k);  
                    cons(j).middle =    trace(ac.Y_k(k)*W_f)       ...
                                        + ac.P_D(t, k)                ...
                                        - ac.C_w(k)*wind.P_wf(t);
                    cons(j).upper =     ac.P_max(k);
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(1) P_inj';
                    
                    % (2) Q_inj
                    j = j+1;
                    cons(j).lower =     ac.Q_min(k);  
                    cons(j).middle =    trace(ac.Ybar_k(k)*W_f)    ...
                                        + ac.Q_D(t, k);
                    cons(j).upper =     ac.Q_max(k);
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(2) Q_inj';
                    
                    % (3) V_bus
                    j = j+1;
                    cons(j).lower =     ac.V_min(k)^2;  
                    cons(j).middle =    trace(ac.M_k(k)*W_f);
                    cons(j).upper =     ac.V_max(k)^2;
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(3) V_bus';
                    
                    
                    
                end
                
                for i = 1:N
                
                    for k = 1:ac.N_b

                        % (1") P_inj scen
                        j = j + 1;
                        cons(j).lower =     ac.P_min(k);  
                        cons(j).middle =    trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                                            + ac.P_D(t, k)                 ...
                                            - ac.C_w(k)*wind.P_w(t, i);
                        cons(j).upper =     ac.P_max(k);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(1") P_inj';
                                           
                        % (2") Q_inj scen
                        j = j+1;
                        cons(j).lower =     ac.Q_min(k);  
                        cons(j).middle =    trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                                            + ac.Q_D(t, k);
                        cons(j).upper =     ac.Q_max(k);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(2") Q_inj';
                        
                        % (3") V_bus scen
                        j = j+1;
                        cons(j).lower =     ac.V_min(k)^2;  
                        cons(j).middle =    trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, i)));
                        cons(j).upper =     ac.V_max(k)^2;
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(3") V_bus';
                        
                        
                                              
                    end
                    
                    for l = 1:ac.N_G
                        
                        % reserve power is always between R_us and R_ds
                        j = j + 1;
                        k = ac.Gens(l);
                        cons(j).lower =     -1*R_ds(l);
                        cons(j).middle =    trace(ac.Y_k(k)*(W_m*wind.P_m(t,i))) ...
                                             - ac.C_w(k)*wind.P_m(t, i);
                        cons(j).upper =     R_us(l);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(7) bound R';
 
                    end
                    
                    % PSD on W_s
                    j = j + 1;
                    cons(j).middle = W_f + W_m * wind.P_m(t,i);
                    cons(j).description = '(6") PSD-ness W_s';
                    cons(j).type = 2;
                    cons(j).scenario = i;

                end
                
                
                
                % PSD constraints
                j = j + 1;
                cons(j).middle = W_f;
                cons(j).description = '(6) PSD-ness W_f';
                cons(j).type = 2;
                cons(j).scenario = 0;
                               
                % reference bus angle constraints
                j = j + 1;
                cons(j).middle = W_f(k_ref, k_ref);
                cons(j).description = '(5) Refbus angle W_f';
                cons(j).type = 3;
                cons(j).scenario = 0;
                
                % refbus angle on W_s
                j = j + 1;                  
                cons(j).middle = W_m(k_ref, k_ref);
                cons(j).description = '(5") Refbus angle W_m';
                cons(j).type = 3;
                cons(j).scenario = i;
                
                % nonnegativity constraints on R limits
                j = j + 1;
                cons(j).middle = R_us;
                cons(j).description = '(9a) NN R_up';
                cons(j).type = 2;
                cons(j).scenario = 0;
                j = j + 1;
                cons(j).middle = R_ds;
                cons(j).description = '(9b) NN R_down';
                cons(j).type = 2;
                cons(j).scenario = 0;
            elseif strcmpi(obj.name, 'P3*') || strcmp(obj.name, 'P3**')
            %% P3* - with only extreme PSD constraints
                
                varargin = varargin{1};
                ac = varargin{1};   
                wind = varargin{2};
                t = varargin{3};
                decision_vars = varargin{4};

                % decision variables
                W_f = decision_vars{1};
                W_m = decision_vars{2};
                R_us = decision_vars{3};
                R_ds = decision_vars{4};

                % determine number of constraints and preallocate
                [~, N] = size(wind.P_w);
                N_c = 3*ac.N_b + N*(ac.N_b*3 + ac.N_G + 2) + 4;
                [cons(1:N_c).bus] = deal(0); % preallocate struct
                [cons(1:N_c).scenario] = deal(0);
                k_ref = ac.N_b + ac.refbus;
                
                j = 0;
                for k = 1:ac.N_b

                    % (1) P_inj
                    j = j + 1;
                    cons(j).lower =     ac.P_min(k);  
                    cons(j).middle =    trace(ac.Y_k(k)*W_f)       ...
                                        + ac.P_D(t, k)                ...
                                        - ac.C_w(k)*wind.P_wf(t);
                    cons(j).upper =     ac.P_max(k);
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(1) P_inj';
                    
                    % (2) Q_inj
                    j = j+1;
                    cons(j).lower =     ac.Q_min(k);  
                    cons(j).middle =    trace(ac.Ybar_k(k)*W_f)    ...
                                        + ac.Q_D(t, k);
                    cons(j).upper =     ac.Q_max(k);
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(2) Q_inj';
                    
                    % (3) V_bus
                    j = j+1;
                    cons(j).lower =     ac.V_min(k)^2;  
                    cons(j).middle =    trace(ac.M_k(k)*W_f);
                    cons(j).upper =     ac.V_max(k)^2;
                    cons(j).bus = k;
                    cons(j).scenario = 0;
                    cons(j).type = 1;
                    cons(j).description = '(3) V_bus';

                end
                
                for i = 1:N
                
                    for k = 1:ac.N_b

                        % (1") P_inj scen
                        j = j + 1;
                        cons(j).lower =     ac.P_min(k);  
                        cons(j).middle =    trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                                            + ac.P_D(t, k)                 ...
                                            - ac.C_w(k)*wind.P_w(t, i);
                        cons(j).upper =     ac.P_max(k);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(1") P_inj';
                                           
                        % (2") Q_inj scen
                        j = j+1;
                        cons(j).lower =     ac.Q_min(k);  
                        cons(j).middle =    trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                                            + ac.Q_D(t, k);
                        cons(j).upper =     ac.Q_max(k);
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(2") Q_inj';
                        
                        % (3") V_bus scen
                        j = j+1;
                        cons(j).lower =     ac.V_min(k)^2;  
                        cons(j).middle =    trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, i)));
                        cons(j).upper =     ac.V_max(k)^2;
                        cons(j).bus = k;
                        cons(j).scenario = i;
                        cons(j).type = 1;
                        cons(j).description = '(3") V_bus';                       
                                              
                    end
                    
                    if strcmp(obj.name, 'P3*')
                        for l = 1:ac.N_G

                            % reserve power is always between R_us and R_ds
                            j = j + 1;
                            k = ac.Gens(l);
                            cons(j).lower =     -1*R_ds(l);
                            cons(j).middle =    trace(ac.Y_k(k)*(W_m*wind.P_m(t,i))) ...
                                                 - ac.C_w(k)*wind.P_m(t, i);
                            cons(j).upper =     R_us(l);
                            cons(j).bus = k;
                            cons(j).scenario = i;
                            cons(j).type = 1;
                            cons(j).description = '(7) bound R';

                        end
                    end
                    
                end
                
                % PSD constraints
                j = j + 1;
                cons(j).middle = W_f;
                cons(j).description = '(6) PSD-ness W_f';
                cons(j).type = 2;
                cons(j).scenario = 0;
                
                % PSD on W_s min
                j = j + 1;
                cons(j).middle = W_f + min(wind.P_m(t, :))*W_m;
                cons(j).description = '(6") PSD-ness W_s min';
                cons(j).type = 2;
                cons(j).scenario = 0;
                
                % PSD on W_s max
                j = j + 1;
                cons(j).middle = W_f + max(wind.P_m(t, :))*W_m;
                cons(j).description = '(6") PSD-ness W_s max';
                cons(j).type = 2;
                cons(j).scenario = 0;

                % reference bus angle constraints
                j = j + 1;
                cons(j).middle = W_f(k_ref, k_ref);
                cons(j).description = '(5) Refbus angle W_f';
                cons(j).type = 3;
                cons(j).scenario = 0;
                
                % refbus angle on W_s
                j = j + 1;                  
                cons(j).middle = W_m(k_ref, k_ref);
                cons(j).description = '(5") Refbus angle W_m';
                cons(j).type = 3;
                cons(j).scenario = i;
                
                % nonnegativity constraints on R limits
                j = j + 1;
                cons(j).middle = R_us;
                cons(j).description = '(9a) NN R_up';
                cons(j).type = 2;
                cons(j).scenario = 0;
                j = j + 1;
                cons(j).middle = R_ds;
                cons(j).description = '(9b) NN R_down';
                cons(j).type = 2;
                cons(j).scenario = 0;
                
            elseif strcmpi(obj.name, 'DC')
            %% DC
            
                varargin = varargin{1};
                dc = varargin{1};   
                wind = varargin{2};
                t = varargin{3};
                decision_vars = varargin{4};

                % decision variables
                P_G = decision_vars{1};
                R_us = decision_vars{2};
                R_ds = decision_vars{3};
                d_us = decision_vars{4};
                d_ds = decision_vars{5};

                % determine number of constraints and preallocate
                [~, N] = size(wind.P_w);
                N_c = 1;
                [cons(1:N_c).bus] = deal(0); % preallocate struct
                [cons(1:N_c).scenario] = deal(0);
                
                j = 0;
           
                cons(1).lower =     dc.P_Gmin;
                cons(1).middle =    P_G;
                cons(1).upper =     dc.P_Gmax;
                cons(1).bus = 0;
                cons(1).scenario = 0;
                cons(1).type = 1;
                cons(1).description = 'Bounds P_G';
                
                
                 % define deterministic power injection vector  
                P_injf = (dc.C_G * P_G - dc.P_D(t, :)' + dc.C_w * wind.P_wf(t));

                % power balance constraints
                cons(2).middle =    ones(1, dc.N_b) * P_injf;
                cons(2).bus = 0;
                cons(2).scenario = 0;
                cons(2).type = 3;
                cons(2).description = 'Power balance';
                
                % line flow limits
                cons(3).lower =     - dc.P_fmax;
                cons(3).middle =    dc.B_f * [dc.B_bustildeinv * P_injf(1:end-1); 0];
                cons(3).upper =     dc.P_fmax;
                cons(3).bus = 0;
                cons(3).scenario = 0;
                cons(3).type = 1;
                cons(3).description = 'Line flow limits';
                
                j = 3;
                % loop over scenarios
                for i = 1:N

                    % define reserve power
                    R = d_us * max(0, -wind.P_m(t, i)) - d_ds * max(0, wind.P_m(t, i));

                    % define scenario power injection vector
                    P_injs = dc.C_G * (P_G + R) + dc.C_w * wind.P_w(t, i) - dc.P_D(t, :)';

                    % power balance constraints
                    j = j+1;
                    cons(j).middle =    ones(1, dc.N_b) * P_injs;
                    cons(j).bus = 0;
                    cons(j).scenario = i;
                    cons(j).type = 3;
                    cons(j).description = 'Power balance';

                    % generator limits
                    j = j + 1;
                    cons(j).lower =     dc.P_Gmin;
                    cons(j).middle =    P_G + R;
                    cons(j).upper =     dc.P_Gmax;
                    cons(j).bus = 0;
                    cons(j).scenario = i;
                    cons(j).type = 1;
                    cons(j).description = 'Bounds P_G+R';

                    % line flow limits
                    j = j + 1;
                    cons(j).lower =     - dc.P_fmax;
                    cons(j).middle =    dc.B_f * [dc.B_bustildeinv * P_injs(1:end-1); 0];
                    cons(j).upper =     dc.P_fmax;
                    cons(j).bus = 0;
                    cons(j).scenario = i;
                    cons(j).type = 1;
                    cons(j).description = 'Line flow limits';
                               
                    % reserve requirements constraints
                    j = j + 1;
                    cons(j).lower =     - R_ds;
                    cons(j).middle =    R;
                    cons(j).upper =     R_us;
                    cons(j).bus = 0;
                    cons(j).scenario = i;
                    cons(j).type = 1;
                    cons(j).description = 'Bounds R';
                    
                   
                end
                
                j = j + 1;
                cons(j).middle = R_us;
                cons(j).description = '(9a) NN R_up';
                cons(j).type = 2;
                cons(j).scenario = 0;
                cons(j).bus = 0;
                
                j = j + 1;
                cons(j).middle = R_ds;
                cons(j).description = '(9b) NN R_down';
                cons(j).type = 2;
                cons(j).scenario = 0;
                cons(j).bus = 0;
                
            end
            % ==============================  
            % Insert other formulations here
            % ==============================
        end
        
        function Obj = get_objective(obj, varargin)
        %% returns the objective function
        %
        % Arguments
        % ---------
        % see prepare() or evaluate() 
            
            Obj = 0;

            if strcmpi(obj.name, 'P1')
            %% P1
                
                varargin = varargin{1};
                ac = varargin{1};   
                wind = varargin{2};
                t = varargin{3};
                decision_vars = varargin{4};
                [~,N] = size(wind.P_w);
                
                
                % decision variables
                W_f = decision_vars{1};
                W_s = decision_vars{2};
                R_us = decision_vars{3};
                R_ds = decision_vars{4};
                
                if length(varargin) < 5
                    lambda = 0.5;
                else
                    lambda = varargin{5};
                end

                
                % P_G for every generator
                for j = 1:ac.N_G
                    k = ac.Gens(j);
                    Obj = Obj + ac.c_us(j)*trace(ac.Y_k(k)*W_f);
                    for i = 1:N
                        Obj = Obj + 1/N * ac.c_us(j) * trace(ac.Y_k(k)*W_s(:,:,i));
                    end
                end
                Obj = Obj + lambda*(ac.c_us' * R_us + ac.c_ds' * R_ds);
                                    
            elseif strcmpi(obj.name, 'P3') || strcmp(obj.name, 'P3*')
                %% P3*
                
                varargin = varargin{1};
                ac = varargin{1};   
                wind = varargin{2};
                t = varargin{3};
                decision_vars = varargin{4};
                
                % decision variables
                W_f = decision_vars{1};
                W_m = decision_vars{2};
                R_us = decision_vars{3};
                R_ds = decision_vars{4};

                if length(varargin) < 5
                    lambda = 0.5;
                else
                    lambda = varargin{5};
                end
                
                % P_G for every generator
                for j = 1:ac.N_G
                    k = ac.Gens(j);
                    Obj = Obj + ac.c_us(j)*trace(ac.Y_k(k)*W_f);
                end
                
                %  Reserve requirements
                Obj = Obj + lambda*(ac.c_us' * R_us + ac.c_ds' * R_ds);
                
            elseif strcmpi(obj.name, 'P3**')
            %% P3** with weird objective
                
                varargin = varargin{1};
                ac = varargin{1};   
                wind = varargin{2};
                t = varargin{3};
                decision_vars = varargin{4};
                
                % decision variables
                W_f = decision_vars{1};
                W_m = decision_vars{2};
                R_us = decision_vars{3};
                R_ds = decision_vars{4};

                if length(varargin) < 5
                    lambda = 0.5;
                else
                    lambda = varargin{5};
                end
                
                % P_G for every generator
                for j = 1:ac.N_G
                    k = ac.Gens(j);
                    Obj = Obj + ac.c_us(j)*trace(ac.Y_k(k)*W_f);
                    
                    % up and downspinning costs
                    mismatch_and_zero = [wind.P_m(t,:) 0];
                    Obj = Obj + lambda * (ac.c_us(j)*trace(ac.Y_k(k)*min(mismatch_and_zero)*W_m) ...
                              + ac.c_us(j)*trace(ac.Y_k(k)*max(mismatch_and_zero)*W_m));  
                end
                
                % Reserve requirements
%                 Obj = Obj + lambda*(ac.c_us' * R_us + ac.c_ds' * R_ds);
                
                
            elseif strcmpi(obj.name, 'DC')
            %% DC
                varargin = varargin{1};
                dc = varargin{1};   
                decision_vars = varargin{4};
                
                % decision variables
                P_G = decision_vars{1};
                R_us = decision_vars{2};
                R_ds = decision_vars{3};
                
                % loop over generators and add generator cost at time t
                for k = 1:dc.N_G
                    Obj = Obj + dc.c_qu(k) * (P_G(k))^2 + ...
                                            dc.c_li(k) * P_G(k);
                %       Obj = Obj + dc.c_us(k) * P_G(k);
                end

                % add reserve requirements costs
                Obj = Obj + (dc.c_us' * R_us + dc.c_ds' * R_ds);
            
                
            % ==============================
            % Insert other formulations here
            % ==============================
 
            end
            
        end
        
        function [Objective, Constraints] = prepare(obj, varargin)
        %% returns the constraint set and objective function
        % 
        % Arguments for 'P1' formulation
        % ----------------------------------
        % ac : ac model used
        % wind : wind model used
        % t : current timestep
        % decision_vars : cell with sdpvars
        % lambda : scalar, weight for reserve 
        

        

            % retrieve structure with all the constraints
            cons_struct = obj.get_constraints(varargin);
            
            % convert to LMI for yalmip
            Constraints = [];
            for c = cons_struct

                if c.type == 1      % box constraint
                    Constraints = [Constraints, (c.lower <= c.middle <= c.upper):sprintf('%s | B%3i | S%3i', c.description, c.bus, c.scenario)];
                elseif c.type == 2  % inequality / psd constraint
                    Constraints = [Constraints, (c.middle >= 0):sprintf('%s | B%3i | S%3i', c.description, c.bus, c.scenario)];
                elseif c.type == 3  % equality constraint
                    Constraints = [Constraints, (c.middle == 0):sprintf('%s | B%3i | S%3i', c.description, c.bus, c.scenario)];
                end
            end
            
            % get the objective and return
            Objective = obj.get_objective(varargin);
            
            obj.Constraints = Constraints;
            obj.prepped_constraints = cons_struct;
            obj.Objective = Objective;
        end
        
        function [failed, failed_scen, failed_det] = evaluate(obj, varargin)
        %% evaluate the solution


            % tolerance for constraint satisfaction
            tol = 1e-5;
                
            % get constraints
            cons = obj.get_constraints(varargin);
            obj.checked_constraints = cons; % save
            n_det = length(find([cons.scenario] == 0));
            n_scen = length(cons) - n_det;

            % evaluate limits with tolerance
            failed_det = 0;
            failed_scen = 0;


            printout = '| BUS | SCEN | DESCRIPTION \t \t \t \t |  LOWER  |  ACTUAL |  UPPER  |\n|======================================================================|\n';

            % check all the constraints
            for c = cons

                if c.type == 1      % box constraint 
                    if any(c.lower - tol >= c.middle) || any(c.middle >= c.upper + tol)
                        if c.scenario == 0
                            failed_det = failed_det + 1;
                        else
                            failed_scen = failed_scen + 1;
                        end
                        printout = [printout sprintf('| %3i | %4i | %25s | %7.4f | %7.4f | %7.4f |\n', ...
                         c.bus, c.scenario, c.description, c.lower, c.middle, c.upper)];
                    end
                elseif c.type == 2  % PSD constraint
                    [~,n] = size(c.middle);
                    if ((n==1) && any(c.middle < -tol))  ...
                    || ((n> 1) && any(eig(c.middle) < -tol))
                        if c.scenario == 0
                            failed_det = failed_det + 1;
                        else
                            failed_scen = failed_scen + 1;
                        end
                        if n > 1
                            c.middle = min(eig(c.middle));
                        else
                            c.middle = min(c.middle);
                        end
                        printout = [printout sprintf('| %3i | %4i | %25s | %7.4f | %7.4f | %7.4f |\n', ...
                         c.bus, c.scenario, c.description, c.lower, c.middle, c.upper)];
                    end

                elseif c.type == 3  % equality constraint
                    if any(abs(c.middle) > tol)
                        if c.scenario == 0
                            failed_det = failed_det + 1;
                        else
                            failed_scen = failed_scen + 1;
                        end
                        printout = [printout sprintf('| %3i | %4i | %25s | \t %15.4d is not 0  |\n', ...
                         c.bus, c.scenario, c.description, c.middle)];
                    end

                end
            end
            printout = [printout '|======================================================================|\n\n'];

            % output
            
            failed = failed_det + failed_scen;
            if failed
                if failed_det
                    fprintf('%i/%i determisitic constraints violated\n', failed_det, n_det);
                end
                if failed_scen
                    fprintf('%i/%i scenario constraints violated\n', failed_scen, n_scen);
                end
                disp('<a href="matlab:p.print_details">Show details</a>');
                obj.details = printout;
            else
                fprintf('All constraints satisfied\n');
            end
            
            
        end
        
        function print_details(obj)
        %% prints out the details of the evaluation to the command line
            assert(~isempty(obj.details), 'Run evaluation first');
            fprintf(['\n\n' obj.details]);
        end
        
    end
    
end
    
    
