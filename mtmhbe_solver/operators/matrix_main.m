function M = matrix_main(xsec, grid, paths)
    
    % Initialize grids
    grid = energy_grid(grid);
    
    % Load Cross Sections Files
    xsec = cross_sections_main(xsec, paths);
    
    % Initialize z-index
    zin = init_z_index(xsec, grid);
    
    %% Generate matrix operator
    
    % Initalize Operator Structure M
    M = initialize_matrix(xsec, grid);
    
    % Coupling Operators
    M = harmonic_operator(M, grid, zin);
    M = ionization_operator(M, grid, zin);
    M = field_operator(M, grid, zin);
    
    % Collision Operators
    [M, xsec, zin] = collision_operators(M, grid, zin, xsec);
    
    % QSS Mass Constraint
    mass.I = grid.INeps(:)*0+1;
    mass.J = grid.INeps(:);
    mass.V = grid.E_INT_MASS(:);
    mass.Z = repmat(zin.mass_int, grid.Neps, 1);
    M = append_IJVZ(M, mass, [], []);
    
    % Reduce I,J,V,Z
    k = logical(M.I);
    M.I = M.I(k);
    M.J = M.J(k);
    M.Z = M.Z(k);
    M.V = M.V(k);
    
    % QSS Cpz
    [IJ, ~, ic] = unique([M.J, M.I], 'rows', 'sorted');
    IJ = cast(IJ, M.int_type);
    ic = cast(ic, M.int_type);
    M.Iu = IJ(:, 2); % unique I in [I, J]
    M.Ju = IJ(:, 1); % unique J in [I, J]
    M.Cpz = sparse(ic, M.Z, M.V, size(IJ, 1), zin.Nz);
    
    % Get inidices of indices to zero in QSS Case
    [ii, jj] = find(M.Cpz);
    kk = find(M.Iu(ii)==1 & jj~=zin.mass_int);
    M.ind_bc_zero = sub2ind(size(M.Cpz), ii(kk), jj(kk));
    
    %% Optional: get statistics on C_pz co-location
    calc_jac_sizes = true;
    M.Nz_act = NaN;
    M.stacked_frac = NaN;
    if calc_jac_sizes
    
        % Temp copy of V
        VV = M.V;
        VV(M.I==1 & M.Z~=zin.mass_int) = 0; % Clear elements that are overwritten by mass conservation
        Nnz_raw = nnz(VV(:));
    
        % Clear bc zero from Y, assuming QSS model
        YY = M.Cpz;
        YY(M.ind_bc_zero) = 0.0;
    
        % Estimate of Raw and Actual Cpz entries
        M.Nz_act = nnz(YY);
        M.stacked_frac = (Nnz_raw-M.Nz_act)./Nnz_raw;
    
        % Example Jac to compare
        ZZZ = sparse(M.Iu, M.Ju, sum(YY, 2), M.N, M.N);
    
        M.jac_entries_k = nnz(ZZZ) / 1000.0;
        M.jac_coll_den  = nnz(ZZZ(1:grid.Neps, 1:grid.Neps))/(double(grid.Neps)^2) * 100; % Percentage
        M.Cpz_entries_k = M.Nz_act / 1000.0;
        M.Cpz_coincident = M.stacked_frac * 100;
    
    end
    
    %% Finalize matrix structure
    M = finalize_matrix(M, grid, zin, xsec);
    
end
    
    
%% Fill in initial z-index data (not related to species)
function zin = init_z_index(xsec, grid)
    
    % Indicator fields, these state where in z each category is located
    zin.mass_int = 1;
    zin.field = 2;
    zin.omega = 3;
    zin.nubar = zin.omega(end) + double(1:grid.NL0)';
    zin.coll = zin.nubar(end) + double(1:xsec.Neps)'; % In all cases, at least store user species
    
    % Index Field, also acts as a logical
    zin.u_z = [zeros(zin.nubar(end), 1); double(1:xsec.Neps)']; % Starts
    zin.Nz = numel(zin.u_z);
    zin.ISM = xsec.ISM;
    
    % Cell array for storing inverse s_z array (negative values indicate z
    % holds superelastic process s)
    zin.s_z = cell(zin.Nz, 1);
    zin.s_z(zin.field) = {'Field'};
    zin.s_z(zin.omega) = {'Omega'};
    for i = 1:grid.NL0
        zin.s_z{zin.nubar(i)} = [grid.K(i), grid.R(i)];
    end
    
    % If Bolsig+ model, initialize empty de, grat
    if xsec.ISM == 1 || xsec.ISM == 0
        zin.de = zeros(zin.Nz, 1); % Starts
        zin.gratio = zeros(zin.Nz, 1); % Starts
        zin.is_upper_state = false(zin.Nz, 1);
        zin.is_thermal = false(zin.Nz, 1);
    end
    
end
    

%% Initialize empty matrix operator
function M = initialize_matrix(xsec, g)

    % Store constants
    M.const = boltz_constants;
    
    % Td estimate comes from assuming dense NL0 and tri-banded coupling
    Nest = g.Neps*g.Neps*g.NL0*1.2 + 3 * g.Neps*g.NL0*(g.NLK - g.NL0);
    
    % Determine integer type needed for size of Ax=B
    if Nest < intmax("int16")
        itype = 'int16';
    elseif Nest < intmax("int32")
        itype = 'int32';
    else
        itype = 'int64';
    end
    % itype = 'double';
    
    % Quasi-Steady State estimate comes from assuming tri-banded on all,
    % block bandwidth of NLKR
    M.I = zeros(Nest, 1, itype);
    M.J = zeros(Nest, 1, itype);
    M.Z = zeros(Nest, 1, itype);
    M.V = zeros(Nest, 1);
    M.N = g.N;
    M.Neps = g.Neps;
    M.row_cnt = 1;
    M.int_type = itype;
    
    % Store empty rate matrix
    M.rates.integral = zeros(xsec.Ns, g.Neps);
    M.rates.zid = zeros(xsec.Ns, 1);

end


%% Finalize matrix
function M = finalize_matrix(M, grid, zin, xsec)

    % Clear original pointers
    M = rmfield(M, 'I');
    M = rmfield(M, 'J');
    M = rmfield(M, 'Z');
    M = rmfield(M, 'V');
    
    % Store in M
    M.grid = grid;
    M.zin = zin;
    M.xsec = xsec;

end
