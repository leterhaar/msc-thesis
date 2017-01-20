function PGvsPW(type, varargin)
    check_class(type, 'char');
    
    %% extract settings from caller
    wind = evalin('caller', 'wind');
    ac = evalin('caller', 'ac');
    t = evalin('caller', 't');
    
    %% define PGs
    PGs = nan(3, ac.N_G);
    Pwf = wind.P_wf(t);
    Pms = [min(wind.P_m(t,:)) 0 max(wind.P_m(t,:))];
    
    if strcmpi(type, 'p3')
        W_0 = varargin{1};
        W_us = varargin{2};
        W_ds = varargin{3};
        
        for i = 1:length(Pms)
            W_s = W_0 + min(Pms(i) + Pwf, Pwf)*W_us + ...
                max(Pms(i), 0)*W_ds;
            
            for j = 1:ac.N_G
                k = ac.Gens(j);
                PGs(i,j) = trace(ac.Y_k(k)*W_s) - ac.C_w(k)*(Pms(i) + Pwf) + ac.P_D(t, k);
            end
        end
    elseif strcmpi(type, 'p4')
        W_max = varargin{1};
        W_us = varargin{2};
        W_ds = varargin{3};
        Pwmax = varargin{4};
        
        for i = 1:length(Pms)
            W_s = W_max + max(-Pms(i), 0) * W_us + ...
                  min(-(Pms(i)+Pwf) + Pwmax, -Pwf + Pwmax) * W_ds;
              
            
            for j = 1:ac.N_G
                k = ac.Gens(j);
                PGs(i,j) = trace(ac.Y_k(k)*W_s) - ac.C_w(k)*(Pms(i) + Pwf) + ac.P_D(t, k);
            end
        end
        
    elseif strcmpi(type, 'p1')
        W_ss = varargin{1};
        W_f = varargin{2};
        
        for i = 1:length(Pms)
            if Pms(i) == 0
                W_s = W_f;
            else
                i_scen = find(wind.P_m(t,:) == Pms(i), 1);
                W_s = W_ss(:,:,i_scen);
            end
            
            for j = 1:ac.N_G
                k = ac.Gens(j);
                PGs(i,j) = trace(ac.Y_k(k)*W_s) - ac.C_w(k)*(Pms(i) + Pwf) + ac.P_D(t, k);
            end
        end
    else
        error('Unknown type: %s', type);
    end
    
    % something is wrong here.... 
    
    % plot
    initfig('P^G + R vs P^m', str2double(type(end)));
    plot(repmat(Pms', 1, ac.N_G), PGs);
    xlim([-Pwf max(Pms)]);
    gennames = strsplit(sprintf('G%i|', ac.Gens), '|');
    legend(gennames{1:end-1});

end
