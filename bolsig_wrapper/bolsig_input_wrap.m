
%% Wrapper function for converting p/g structs to a common set of input structs for Bolsig+
function bdat = bolsig_input_wrap(xsec, paths, gas, field, grid, settings)

    % READCOLLISIONS, in coll struct 
    const = boltz_constants;
    xsec = cross_sections_main(xsec, paths);
    bdat.file_names =  xsec.file_paths;
    % bdat.file_names =  xsec.file_names;

    % use all species for each file, bolsig plays nice with this
    bdat.species_names = xsec.spec_names; 
    bdat.extrapolate = xsec.extrap;
    
    % Conditions
    bdat.eb_cosine_field_angle = 0.0; % zero for any AC field
    bdat.gas_temp_K = gas.Tgas; % zero for any AC field
    bdat.excitation_temp_K = gas.Texc; % zero for any AC field
    bdat.transition_energy = 0.0; % When zero, super-elastic collision computed by boltzmann distribution
    bdat.ionization_degree = gas.ne_N; 
    bdat.gas_pressure = gas.press_Pa;
    if bdat.gas_temp_K==0
        bdat.N0 = bdat.gas_pressure ./ (const.KB .* 300.0);
    else
        bdat.N0 = bdat.gas_pressure ./ (const.KB .* bdat.gas_temp_K);
    end
    bdat.reduced_efield_rms = field.EN_TD; % convert amplitude
    if isfield(field, 'omega') && field.omega~=0.0
        bdat.reduced_angular_frequency = field.omega ./ bdat.N0; % For AC field=
        bdat.reduced_efield_rms = bdat.reduced_efield_rms / sqrt(2.0);
    else
        bdat.reduced_angular_frequency = 0.0; % For DC field
    end
    bdat.electron_density = bdat.N0 * gas.ne_N; % Electron number density for Coulomb logarithm
    bdat.ion_charge_param = 0; % Ignores ion effects, valid out to ne/N of about ne/N<1e-4 for E/N=10Td, and ne/N<1e-3 for E/N<100Td
    bdat.ion_neutral_mratio = 1; % not used if ion effects ignored
    bdat.coulomb_effects = 0; % default to 0, No F1 term or modified log.
    if settings.coulomb_collision_mod
        bdat.coulomb_effects = 2; % Add modified logarithm, no F1
    end 
    if settings.coulomb_f1_component % add on F1 
        bdat.coulomb_effects = bdat.coulomb_effects + 1;
    end
    bdat.ionization_sharing = 1; % share ionized electron energy (1=equal, 2= one takes all energy)
    bdat.problem_type = 1; % Typical run type for bolsig+ calculations
    bdat.mawellian_mean_energy = 0; % Does not assume electron temperature

    NN = grid.Neps;
    if settings.bolsig_quick_run
        NN = min(1000, grid.Neps);
    end
    bdat.grid_count = NN;
%     bdat.grid_count = grid.Neps;
    bdat.manual_grid = settings.bolsig_grid_type; % log?
    bdat.max_energy_eV = grid.eV_max;
    bdat.precision = settings.tol_abs;
    toll = settings.tol_rel;
    if settings.bolsig_quick_run
        toll = max(1e-6, settings.tol_rel);
        bdat.precision = max(1e-20, settings.tol_abs);
    end
    bdat.convergence = toll;
    bdat.max_iterations = 2000;
    if iscell(gas.spec_frac) % Store cell array of all species fractions
        bdat.species_fractions = gas.spec_frac;
    else % use old method, for stability
        bdat.species_fractions{1} = gas.spec_frac;
%         Nspecs = numel(gas.spec_frac);
%         bdat.species_fractions = cell(Nspecs, 1);
%         for i = 1:Nspecs
%             bdat.species_fractions{i} = gas.spec_frac(i);
%         end
    end

    bdat.normalize_species_fractions = 1;

end
