function M = field_operator(M, g, zin)
    
    const = boltz_constants;

    % Compute L1 and L2 Coefficients
    [L1, L2] = L1_L2_coeff(g.L, g.L');
    
    % Delta coefficient
    Delta_KpK = Delta_KpK_coeff(g.K, g.R, g.K', g.R');
    
    % Merge Delta with L1, sets non-zero entries
    Delta_k_L1 = Delta_KpK .* L1;
    
    % Find Non-Zero in L1/L2, Delta
    ind = find( logical(Delta_k_L1)  );
    [i, j] = ind2sub(size(L1), ind);
    
    % Reshape Coefficients to dim 3
    Delta_k_L1 = reshape(Delta_k_L1(ind), 1, 1, []);
    L2 = reshape(L2(ind), 1, 1, []);
    
    % Boolean for even entries along dim 3 (matching coeff)
    L_is_even = double(~mod(g.L(i), 2));
    L_is_even = reshape(L_is_even, 1, 1, []);
    
    
    %% Even grid, effect from odd
    eo_m0 = g.OC;
    eo_m1 = [0; g.OC(1:end-1)];
    eo_m2 = [0; 0; g.OC(1:end-2)];
    Ap = (0.5.*(g.ER-g.EL).*(L2 - 1.0) - g.EL)./(eo_m1 - eo_m2);
    Bp = (0.5.*(g.ER-g.EL).*(L2 - 1.0) + g.ER)./(eo_m0 - eo_m1);
    Fm2 = Ap.*(eo_m1-g.EL);
    Fm1 = Ap.*(g.EL-eo_m2) + Bp.*(eo_m0-g.ER);
    Fm0 = Bp.*(g.ER-eo_m1);
    
    Fm2(1) = 0;
    Fm1(1) = Bp(1).*(eo_m0(1)-g.ER(1));
    Fm0(1) = Bp(1).*(g.ER(1)-eo_m1(1));
    
    % Cat to Neps x 3
    Ve = [Fm2, Fm1, Fm0] ./ g.Evol_0;
    Je = [g.INeps-2,  g.INeps-1, g.INeps];
    
    % Apply BC's
    Je(1, 1:2, :) = 1;   % Make index valid
    Ve(1, 1:2, :) = 0.0; % Clear entry
    Je(2,   1, :) = 1;   % Make index valid
    Ve(2,   1, :) = 0.0; % Clear entry
    Je = repmat(Je, 1, 1, numel(ind));
    

    %% Generate odd grid from anisotropy on even grid
    ee_m1 = [0; g.EC(1:end-1)];
    ee_m0 = g.EC;
    ee_p1 = [g.EC(2:end); 2*g.EC(end)-g.EC(end-1)];
    Bp = ( 0.5.*(g.OR-g.OL).*(L2 - 1.0) + g.OR)./(ee_p1 - ee_m0);
    Ap = ( 0.5.*(g.OR-g.OL).*(L2 - 1.0) - g.OL)./(ee_m0 - ee_m1);
    
    Fm1 = Ap.*(ee_m0-g.OL);
    Fm0 = Ap.*(g.OL-ee_m1) + Bp.*(ee_p1-g.OR);
    Fp1 = Bp.*(g.OR-ee_m0);
    
    % Cat to Neps x 3
    Vo = [Fm1, Fm0, Fp1] ./ g.Ovol_0;
    Jo = [g.INeps-1, g.INeps, g.INeps+1];
    
    % High Energy BC on odd L, references even grid at eps_max + 1/2
    Jo(1, 1, :) = Jo(1, 2, :); % Assume even solution is constant at left edge, very small effect
    Jo(g.Neps, 3, :) = Jo(g.Neps, 2, :) ; % Make index valid
    Jo = repmat(Jo, 1, 1, numel(ind));
    
    %% Combine odd and even
    D.I = repmat(g.INeps, 1, 3, numel(ind));
    D.J = Je.*L_is_even + Jo.*(1.0-L_is_even);      
    D.Z = D.I.*0 + zin.field;
    D.V = Ve.*L_is_even + Vo.*(1.0-L_is_even);
    D.V = - const.GAMMA .* Delta_k_L1 .* D.V;
    
    % Shift dim3 to correct i/j tile
    i_tile = reshape((i-1).*g.Neps, 1, 1, []);
    D.I = D.I + i_tile;
    j_tile = reshape((j-1).*g.Neps, 1, 1, []);
    D.J = D.J + j_tile;
    
    % Append to Stack
    if ~isempty(D)
        M = append_IJVZ(M, D, [], []);
    end

end
