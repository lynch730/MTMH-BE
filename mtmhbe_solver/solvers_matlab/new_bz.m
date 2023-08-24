function [b_z, nubar, wtime] = new_bz(t, X, p, M)
    
    % Start Timer
    bz_time = tic;
    
    % Number Density
    Tgas = p.gas.Tgas;
    Tref = Tgas;
    if Tref<=0
        Tref = 300.0;
    end
    Ngas = p.gas.press_Pa ./ (M.const.KB .* Tref);
    
    % Species Densities, assumed to be independent of time here
    b_z = p.gas.N_z;
    
    % Thermal Term
    KbTe = Tgas .* M.const.KB ./ M.const.QE;
    b_z(M.zin.is_thermal) = b_z(M.zin.is_thermal) .* KbTe;
    
    % Get Efield either from function or constant
    if isa(p.EN_TD,'function_handle')
        EN_TD = p.EN_TD(t);
    else
        EN_TD = p.EN_TD;
    end
    
    % Convert E/N in Td to E in V/m
    b_z(M.zin.field) = EN_TD .* M.const.VMM_PER_TD .* Ngas;
    
    % Attachment
    nubar = zeros(M.grid.NL0, 1);
    if M.zin.nubar_eval.any_nonc
        
         % Reshape X to L=0
        FL0 = reshape(X(1:M.grid.Neps*M.grid.NL0), M.grid.Neps, M.grid.NL0);
        
        % Compute Densities
        N_z_tmp = reshape(p.gas.N_z(M.zin.nubar_eval.z_s), [], 1);
        N_z_tmp(M.zin.nubar_eval.ind_att) = -N_z_tmp(M.zin.nubar_eval.ind_att);
        
        % Compute Total Nepsbar
        nubar = reshape(sum(N_z_tmp .* M.zin.nubar_eval.integral, 1) * FL0, [], 1);
        
    end
    
    % Apply b_z
    b_z(M.zin.nubar) = nubar;
    
    % Field Frequencies
    if M.grid.NK > 1
        b_z(M.zin.omega) = p.omega;
    end
    
    % If time not given, then apply mass constraint (qss)
    if isempty(t)
        b_z(M.zin.mass_int) = 1.0;
    end
    
    % End Timer
    wtime = toc(bz_time);
    
end
   