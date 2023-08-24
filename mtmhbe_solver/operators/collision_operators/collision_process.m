function [sink, source, u] = collision_process(proc, u)
    
    % Get process-specific-grids
    u = process_specific_grids(proc, u);
    
    % defaults
    sink = [];
    source = [];
    
    % Load constants
    const = boltz_constants;
    
    %% Generic function for integrating non-constant cell values
    % specifically for collision terms with dependence on sqrt(eps) in nu
    
    % Local copy of tabulated energy locations
    eps_j = proc.eps_j;
    N_eps_j = numel(eps_j);
    
    % Create grid eps_k
    eps_k = [eps_j; u.eps_i_shift; u.eps_i];
    [eps_k, sort_k, sort_reverse] = unique(eps_k, 'last'); % Last ensures ia points to all OR
    
    % Then we need an index array which states which node each
    % sub-intervals belongs to. This is done by using ia starting at
    % where OL array starts. Then all energy values tied to
    % cross-section's are negative, and made zero. cumsum with logical 
    % makes those indices match the previous
    ind = sort_k - numel(eps_j) - u.N_eps_i_shift;
    ind(ind<1) = 0.0;
    ind = cumsum(logical(ind));
    ind(ind==0) = u.Neps+1;
    ind(end) = [];
    
    %% Cross-Section Solution 
    
    % Main Cross-Section Interpolation
    sigma_k = extended_linear_interp(proc, eps_k);
    
    % Ensure no values less than transition energy
    if proc.is_superelastic
        sigma_k(eps_k<0) = 0.0;
    elseif proc.is_excitation || proc.is_ionization 
        sigma_k(eps_k<=proc.de) = 0.0;
    end
    
    % Recover sigma for the original eps_i data
    u.sigma_i = sigma_k(sort_reverse);
    u.sigma_i = u.sigma_i((N_eps_j + u.N_eps_i_shift + 1):end);
    
    % Store reduced collision frequency
    u.nu_red = const.GAMMA .* sqrt(u.eps_i) .* u.sigma_i;
    
    %% Solve for I_1k, I_2k Sub-Integrals 
    
    % Left face 
    kL = 1:numel(eps_k)-1;
    kR = kL+1;
    
    % Support
    A1 = (eps_k(kR).^1 - eps_k(kL).^1)./1.0; 
    A2 = (eps_k(kR).^2 - eps_k(kL).^2)./2.0; 
    A3 = (eps_k(kR).^3 - eps_k(kL).^3)./3.0;
    A4 = (eps_k(kR).^4 - eps_k(kL).^4)./4.0;
    
    % Cross Section Terms
    M_sig =  (           sigma_k(kR) -            sigma_k(kL)) ./ A1;
    B_sig = -(eps_k(kL).*sigma_k(kR) - eps_k(kR).*sigma_k(kL)) ./ A1;

    % I_k terms
    I_1k = M_sig .* A3 + B_sig .* A2;
    I_2k = M_sig .* A4 + B_sig .* A3;

    %% Accumulate Sink Term Sub Integrals

    % k2i is index array of length of N_eps_k, which stores the ith cell
    % where that k-cell sub-integral belongs to, specifically for sink terms
        % k2i is obtained by exploiting properties of sort_k, whose values
        % are indices of the unsorted eps_k, at positions in the sorted eps_k.
        % As unsorted eps_k has a specific order: eps_j, eps_i_shift, eps_i 
        % subtracting off N_eps_j and eps_i_shift ensures all values <1 
        % are positions in the sorted eps_k containing eps_i. In other
        % words, first two lines below + logical() create a boolean mask
        % for eps_k values that were inserted as eps_i. Those are the
        % fenceposts for sub-integration in the sorted eps_k.
        % Then a new k2i is formed by cumsum of that boolean, where each
        % 1 value increments the following 1,0,0,0,... as indicating a news
        % cell in i to sum to. However, have to be careful with edge cases
        % and off-by-one issues. Also have to create a null bin to throw
        % OOB eps_k values, which is why u.Neps+1 is used and discarded.
    % k2i = sort_k - N_eps_j - u.N_eps_i_shift; % Make all non-eps_i values <=0
    % k2i(k2i<1) = 0.0; % Clear those indices
    % k2i = cumsum(logical(k2i)); % Sum the boolean, creating index array
    % k2i(k2i==0) = u.Neps+1; % Assign OOB values to null
    % k2i(end) = []; % Clear the last value, which is redundant

    % Create I1 by accumulating via k2i array
    I1 = accumarray(ind, I_1k, [u.Neps+1, 1]);
    I1 = I1(1:u.Neps); 
    
    % Create I2 (as I3 in final form)
    I2 = accumarray(ind, I_2k, [u.Neps+1, 1]);
    I2 = I2(1:u.Neps); 
    
    % Assemble I3
    I3 = (I2 - u.eps_c_i.*I1) ./ u.eps_c_i_diff2;
    
    %% Assemble I/J/V for Sink Terms
    
    % Cat matrix with modifiers
    sink.V = - [-I3, I1, I3] .* const.GAMMA;
    
    % Construct Raw Pointers
    sink.I = repmat(u.INeps, 1, 3);
    sink.J = [u.INeps-1, u.INeps, u.INeps+1];
    sink.Z = sink.I*0 + proc.z_s;
    
    % Left BC's
    sink.J(1, 1) = sink.J(1, 2);
    sink.V(1, [1, 3]) = 0.0; 
    
    % Right BC's
    sink.J(u.Neps, 3) = sink.J(u.Neps, 1); % Set to left index
    sink.V(u.Neps, 3) = sink.V(u.Neps, 3) .* u.bc_end_weight;   
    
    % Rate Integral
    sink.integral(1, :) = accumarray( sink.J(:), sink.V(:), [u.Neps, 1] );
    
    % Divide by volume
    sink.V = sink.V ./ u.cell_volume(sink.I);
        
    %% Accumulate Source integrals
    if u.N_eps_i_shift > 0
        
        ind_coll = ind;

        % Row Index - position in sorted energies where 
        % US (UL+de) reside. Zeroed below "1" are the sink energies
        % which are not valid destinations themselves. 
        ind = sort_k - N_eps_j; % Make eps_j and eps_i nodes negative
        [~, ilast] = max(ind); % Sort index of last value for edge case
        ind(ind<1 | ind>u.N_eps_i_shift) = 0; % Preserve only shift values
        ind(ilast) = u.N_eps_i_shift + 1;
        
        % Sum to locate Source bins (j's)
        ind = cumsum(logical(ind));
        ind(ind==0) = u.N_eps_i_shift + 1;
        ind(end) = [];
        
        % Clear OOB
        ind_oob = ind > u.Neps | ind>u.N_eps_i_shift | ind_coll > u.Neps;
        I_1k(ind_oob) = [];
        I_2k(ind_oob) = [];
        ind(ind_oob) = [];
        ind_coll(ind_oob) = [];
        
        % ind and k2i are to be understood as index arrays in 
        % eps_ref aka sub-integral space that indicate the associated
        % source cell (ind) and sink cell (k2i). This means that 
        % in inelastic processes, the sink cell is along the diagonal,
        % and k2i (formally ind) is the row of the sink cell
        % (source.I) with +/- 1 offsets of source.I for columns. 
        % For inelastic sources, that sink is defined at a given J
        % but now applies to a new source row (the energy of the
        % receiving cell). Hence, we use k2i for J and ind for I, 
        % asn ind is now the cell where the sink at J is a source
        %
        % For super-elastics, the sink is associated with a loss to a 
        % higher energy state, so for a given row sink (col), that volume is 
        % a source at a lower energy. Hence the J location of that sink
        % is the same, but the row ind is greater
        
        % Sort into unique row-column locations
        % First column of C is cell i where operator applies
        % Second column of C is cell j where source occurs from
        % ic is the sort array indicating which cell to sum to
        [C, ~, ic] = unique([ind(:), ind_coll(:)], 'rows');
        
        % Sum eps_k
        I1 = accumarray(ic, I_1k, [max(ic), 1]);
        I2 = accumarray(ic, I_2k, [max(ic), 1]);
        
        % Get correct I cell values for I3
        jj = C(:, 2);
        I3 = (I2 - u.eps_c_i(jj).*I1)./ u.eps_c_i_diff2(jj);
        
        % Create Initial I/J/V
        source.I = repmat(C(:, 1), 1, 3);
        source.J = [jj-1, jj, jj+1];
        source.V = [-I3, I1, I3] .* const.GAMMA;
        source.Z = source.I*0 + proc.z_s;
        
        % Left BC
        [ip, ~] = find(source.J<1);
        source.J(ip, 1) = source.J(ip, 2);
        source.V(ip, 1) = 0.0; 
        source.V(ip, 3) = 0.0; 
        
        % Right BC
        [ip, ~] = find(source.J>u.Neps);
        source.J(ip, 3) = source.J(ip, 1); % Set to left index
        source.V(ip, 3) = source.V(ip, 3) .* u.bc_end_weight;   
        
        % Multiply by 2 for ionization
        source.V = source.V .* (1.0 + double(proc.is_ionization));
        
        % Divide by volume
        source.V = source.V ./ u.cell_volume(source.I);

    end
    
end


% Prep odd and even grids for this specific process
%   ue and uo are structs for each odd/even case
function u = process_specific_grids(proc, u)
   
    % Grid Cells
    % u.eps_i = [u.eps_l_i; u.eps_r_i(end)];
    
    % Include Source Energy on even grids (but wont be used for L=0)
    u.eps_i_shift = [];
    if ((proc.is_excitation || proc.is_ionization) && u.is_even)
        u.eps_i_shift = u.eps_i .* (1.0 + double(proc.is_ionization));
        if proc.is_superelastic
            u.eps_i_shift = u.eps_i_shift - proc.de;
        else
            u.eps_i_shift = u.eps_i_shift + proc.de;
        end
        u.eps_i_shift(u.eps_i_shift > u.eps_i(end)) = [];
    end
    
    % Store number of shifted cells
    u.N_eps_i_shift = numel(u.eps_i_shift);

end
