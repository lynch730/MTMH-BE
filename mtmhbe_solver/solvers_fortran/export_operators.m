%% Goal of this function is to export core variabels to mat file
function export_operators(M, paths, fname)
    
    % Rates
    rates_integral = M.rates.integral;
    rates_zid = int32(M.rates.zid(:));
    
    % Grid
    grid_EC = M.grid.EC(:);
    grid_OC = M.grid.OC(:);
    grid_E_INT_ENERGY = M.grid.E_INT_ENERGY(:);
    grid_E_INT_MASS = M.grid.E_INT_MASS(:);
    grid_L_is_0 = int32(M.grid.L_is_0(:));
    
    % Const
    const_KB = M.const.KB;
    const_VMM_PER_TD = M.const.VMM_PER_TD;
    const_QE = M.const.QE;
    
    % Integrals
    integral_attachment = 0.0;
    if isfield(M.integral, 'attachment') > 0
        integral_attachment = M.integral.attachment;
    end
%     integral_attachment = integral_attachment(:);

    integral_ionization = 0.0;
    if isfield(M.integral, 'ionization') > 0
        integral_ionization = M.integral.ionization;
    end
%     integral_ionization = integral_ionization(:);

    % Z index
    zin_N             = int32(M.zin.zid.N);
    zin_nubar         = int32(M.zin.zid.nubar(:));
    zin_elastic_th    = int32(M.zin.zid.elastic_th(:));
    zin_inelastic     = int32(M.zin.zid.inelastic(:));
    zin_superelastic  = int32(M.zin.zid.superelastic(:));
    zin_attach        = int32(M.zin.zid.attach(:));
    zin_ionization    = int32(M.zin.zid.ionization(:));
    zin_field         = int32(M.zin.zid.field(:));
    zin_omega         = int32(M.zin.zid.omega(:));
    zin_new_sid       = int32(M.zin.species.new_sid(:));
    zin_low_type      = int32(M.zin.species.low_type(:));
    zin_ensemble_type = int32(M.zin.species.ensemble_type(:));
    zin_base_specid   = int32(M.zin.species.base_specid(:));
    zin_Zmap          = int32(M.zin.species.Zmap(:));
    zin_mass_int      = int32(M.zin.zid.mass_int(:));
    zin_grat          = M.zin.species.grat(:);
    zin_de            = M.zin.species.de(:);

    % % Cross Section Copy
    xsec_Nproc = int32(M.xsec.Nproc);
    xsec_is_rev = int32(find(M.xsec.is_rev(:)));
    xsec_is_attachment = int32(find(M.xsec.is_attachment(:)));
    xsec_is_excitation = int32(find(M.xsec.is_excitation(:)));
    xsec_is_elastic = int32(find(M.xsec.is_elastic(:)));
    xsec_is_ionization = int32(find(M.xsec.is_ionization(:)));

    % Copy Y Matrix
%     jac_Iu = int32(M.m(:));
%     jac_Ju = int32(M.n(:)); 
    AA = sparse(M.m, M.n, ones(size(M.m)), M.N, M.N);
    [jac_Iu, jac_Ju, ~, sid] = convert_sparse_format(AA);
    
    % if QSS, 
    YY = M.Y(sid, :);
    Iu = M.m(sid);
    [ii, jj] = find(YY);
    kk = find(Iu(ii)==1 & jj~=M.zin.zid.mass_int);
    ind_bc_zero = sub2ind(size(YY), ii(kk), jj(kk));
    [Y_I, Y_J, Y_V, ~, ind_bc_zero] = convert_sparse_format(YY, ind_bc_zero);
    ind_bc_zero = int32(ind_bc_zero);

    % Copy size ints
    N = int32(M.grid.N);
    Neps = int32(M.grid.Neps);
    NL0 = int32(M.grid.NL0);
    
    % M.Y(M.ind_bc_zero) = 0.0;
    % beta = ones(size(M.Y, 2), 1);
    % A1 = sparse(M.m, M.n, M.Y*beta, N, N);
    % A2 = decomposition(A1);
    % A = sparse(M.m, M.n, M.Y*beta*0, N, N);
    % rhs = A * ones(M.N, 1);
    % rhs(1) = 1.0;
    % b = A2\rhs;

    %% Save Integers
    save(fullfile(paths.fortran_solver, fname), ...
        'N', ...                	% I
        'Neps', ...               	% I
        'NL0', ...              	% I
        'rates_integral', ...       % 
        'rates_zid', ...            % 
        'grid_EC', ...              % 
        'grid_OC', ...				% 
        'grid_E_INT_ENERGY', ...    % 
        'grid_E_INT_MASS', ...      % 
        'grid_L_is_0', ...			% 
        'const_KB', ...             % 
        'const_VMM_PER_TD', ...     % 
        'const_QE', ...             % 
        'integral_attachment', ...  % 
        'integral_ionization', ...  % 
        'zin_N', ...                % I
        'zin_nubar', ...			% 
        'zin_elastic_th', ...       % 
        'zin_attach', ...           % 
        'zin_ionization', ...       % 
        'zin_field', ...            % 
        'zin_omega', ...            % 
        'zin_new_sid', ...          % 
        'zin_low_type', ...         % 
        'zin_ensemble_type', ...    % 
        'zin_base_specid', ...      % 
        'zin_Zmap', ...             % 
        'zin_grat', ...             % 
        'zin_de', ...               % 
        'zin_mass_int', ...         %
        'zin_inelastic', ...        %
        'zin_superelastic', ...     %
        'xsec_Nproc', ...           % 
        'xsec_is_rev', ...          % 
        'xsec_is_attachment', ...   % 
        'xsec_is_excitation', ...   % 
        'xsec_is_elastic', ...		% 
        'xsec_is_ionization', ...	% 
        'jac_Iu', ...               %
        'jac_Ju', ...               %
        'Y_I', ...				    % 
        'Y_J', ...				    % 
        'Y_V', ...				    % 
        'ind_bc_zero', ...          %
        '-v7.3');

end



