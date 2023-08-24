function [mdata, mpath] = run_multibolt(paths, M, gas, field, settings, run_name)
    
    EN_TD = field.EN_TD; % convert amplitude to RMS
    if isfield(field, 'omega')
        if field.omega ~= 0.0
            EN_TD = EN_TD / sqrt(2.0);
        end
    end
    
    % xsec = cross_sections_main(xsec, paths);
    % zin = z_index(xsec, 1, grid.NK);
    % [Nz, N0] = species_fractions(gas, M.zin, []);
    % Xz = Nz ./ N0;

    multibolt_folder = paths.multibolt_data;

    % Create file
    if ~exist('run_name','var')
        run_name = datestr(now, 'dd-mm_HH-MM-SS');
    end
    if isunix
        filename = ['multibolt_', run_name, '.sh'];
    else
        filename = ['multibolt_', run_name, '.bat'];
    end
    runfile = fullfile(multibolt_folder, filename);
    fid = fopen( runfile, 'w' );
    
    % Print command
    if isunix
        fprintf(fid, '%s \\', paths.multibolt_exe);
    else
        fprintf(fid, '%s ^', paths.multibolt_exe);
    end
    warg(fid, 'model', 'HD', true);
    warg(fid, 'export_name', run_name);
    warg(fid, 'N_terms', M.grid.NL);
    warg(fid, 'p_Torr', gas.press_Pa .* 0.00750062);
    warg(fid, 'T_K', gas.Tgas);
    warg(fid, 'EN_Td', EN_TD);
    NN = M.grid.Neps;
    if settings.bolsig_quick_run
        NN = min(1000, M.grid.Neps);
    end
    warg(fid, 'Nu', NN);
%     warg(fid, 'Nu', p.grid.Nu);

    warg(fid, 'initial_eV_max', M.grid.eV_max);
    if settings.bolsig_grid_type == 0
        warg(fid, 'USE_ENERGY_REMAP', '');
    end
%     warg(fid, 'use_eV_max_guess', '');
    toll = settings.tol_rel;
    if settings.bolsig_quick_run
        toll = max(1e-6, settings.tol_rel);
    end
    warg(fid, 'conv_err', toll);
    warg(fid, 'weight_f0', 0.8);
    warg(fid, 'iter_max', settings.max_jac_iter);
    warg(fid, 'iter_min', 0);
    if settings.bolsig_grid_type == 0
        warg(fid, 'remap_target_order_span', 10);
        warg(fid, 'remap_grid_trial_max', 20);
    end
    warg(fid, 'multibolt_num_threads', 'max', true);
    warg(fid, 'export_location', multibolt_folder);
    for i = 1:M.xsec.Nfiles
        warg(fid, 'LXCat_Xsec_fid', M.xsec.file_names{i});
    end
%     min_frac = 1e-16;
%     pfrac = gas.spec_frac;
%     pfrac(isnan(pfrac) | pfrac<min_frac) = 0;
%     [~, ind] = max(pfrac);
%     pfrac(ind) = pfrac(ind) + 1.0 - sum(pfrac);
    

    for i = 1:numel(M.xsec.spec_names)
        if gas.spec_frac(i)>0
            wspecies(fid, M.xsec.spec_names{i}, gas.spec_frac(i));
        end
    end


%     wspecies(fid, xsec.pname{1}, 1);
%     wspecies(fid, xsec.pname{2}, 0);
    warg(fid, 'elastic_scattering', 'Isotropic', true);
    warg(fid, 'excitation_scattering', 'Isotropic', true);
    warg(fid, 'ionization_scattering', 'Isotropic', true);
    warg(fid, 'superelastic_scattering', 'Isotropic', true);
    warg(fid, 'sharing', 0.5);
    warg(fid, 'sweep_option', 'EN_Td', true);
    str = sprintf('%10.4f ',  EN_TD);
    wline(fid, ['--sweep_style def ', str]);
    wline(fid, '--DONT_ENFORCE_SUM ');

    % Close fule
    fclose(fid);
    if isunix
        fileattrib(runfile, '+x');
    end

%     com = fileread( fullfile(multibolt_folder, filename) );
%     com = erase(com, '`');
%     com = regexprep(com,'[\n\r]+','');
%     com = [com, ' --help'];
%     system(com);
    % system(sprintf('%s', runfile));
    [mdata.status, mdata.cmdout] = system(sprintf('%s', runfile));

    % get folder path to read data file
    mpath = fullfile(multibolt_folder, run_name);
    
    mdata = extract_multibolt(mpath,  EN_TD);

end

function wspecies(fid, sname, xfrac)
    if xfrac > 1e-4
        x = ['--species "', sname, '" ', num2str(xfrac, '%0.5f')];
    else
        x = ['--species "', sname, '" ', num2str(xfrac, '%4.2e')];
    end
    wline(fid, x)
end

function warg(fid, arg, x, suppress_quotes)
    
    str = ['--', arg, ' '];
    if isnumeric(x)
        if floor(x)==x % isint
            fmt = '%i';
        else
            fmt = '%16.8e';
        end
        x = num2str(x, fmt);
    elseif ~isempty(x) 
        if ~exist('suppress_quotes', 'var')
            suppress_quotes = false;
        end
        if ~suppress_quotes
            x = ['"', x, '"'];
        end
    end
    str = [str, x];
    wline(fid, str)
    
end

function wline(fid, str)
    if isunix 
        fprintf(fid, '\n %s \\', str);
    else
        fprintf(fid, '\n %s ^', str);
    end
end
