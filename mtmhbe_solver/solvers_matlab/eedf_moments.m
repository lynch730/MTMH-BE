function sol = eedf_moments(X, time, M, b_z)
    
    sol.X = X;
    sol.b_z = b_z;
    g = M.grid;
    
    % Time, if used
    Nt = 1;
    if ~isempty(time)
        Nt = time.Nt;
    end
    
    % Load in EEDF
    sol.F = reshape(X, g.Neps, g.NLK, Nt);
    
    %% Mass-Momentum-Energy Moments of EEDF

    % Reshaped F0, F1
    FAe = reshape(reshape(sol.F(:, g.L_is_even, :), g.Neps, [], Nt), g.Neps, []);
    FAo = reshape(reshape(sol.F(:,~g.L_is_even, :), g.Neps, [], Nt), g.Neps, []);
    F1 = reshape(reshape(sol.F(:, g.L==1, :), g.Neps, g.NL1, Nt), g.Neps, []);
    F0 = reshape(reshape(sol.F(:, g.L==0, :), g.Neps, g.NL0, Nt), g.Neps, []);
    
    % Mass-Momentum-Energy
    sol.Fmom.mass( g.L_is_even, :) = reshape(M.grid.E_INT_MASS * FAe, [], Nt);
    sol.Fmom.mass(~g.L_is_even, :) = reshape(M.grid.O_INT_MASS * FAo, [], Nt);
    sol.Fmom.momentum = reshape(M.grid.O_INT_MOM * F1, [g.NL1, Nt]);
    sol.Fmom.energy   = reshape(M.grid.E_INT_ENERGY * F0, [g.NL0, Nt]);
    
    %% Reaction Rates

    % All process rates
    F0 = reshape(F0, g.Neps, [] );
    sol.rates.red.all = reshape( M.rates.integral * F0, [M.xsec.Ns, g.NL0, Nt] );
    
    % Scale by number densities
    Ns = reshape(b_z(M.rates.zid, :), [M.xsec.Ns, 1, Nt]);
    sol.rates.raw.all = Ns .* sol.rates.red.all;
    
    % Raw by type
    sol.rates.raw.superelastic = sol.rates.raw.all(M.xsec.is_superelastic,:,:);
    sol.rates.raw.attachment   = sol.rates.raw.all(M.xsec.is_attachment,:,:);
    sol.rates.raw.excitation   = sol.rates.raw.all(M.xsec.is_excitation,:,:);
    sol.rates.raw.elastic      = sol.rates.raw.all(M.xsec.is_elastic,:,:);
    sol.rates.raw.ionization   = sol.rates.raw.all(M.xsec.is_ionization,:,:);
    
    % Reduced by type
    sol.rates.red.superelastic = sol.rates.red.all(M.xsec.is_superelastic,:,:);
    sol.rates.red.attachment   = sol.rates.red.all(M.xsec.is_attachment,:,:);
    sol.rates.red.excitation   = sol.rates.red.all(M.xsec.is_excitation,:,:);
    sol.rates.red.elastic      = sol.rates.red.all(M.xsec.is_elastic,:,:);
    sol.rates.red.ionization   = sol.rates.red.all(M.xsec.is_ionization,:,:);
    
    % Add grid
    sol.EC = M.grid.EC;
    
end
