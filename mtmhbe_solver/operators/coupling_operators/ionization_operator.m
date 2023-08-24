function M = ionization_operator(M, g, zin)
    
    % Second dimension, LKR columns
    [K, Kp, K1] = ndgrid(g.K, g.K', reshape(g.K(g.L==0), 1, []));
    [R, Rp, R1] = ndgrid(g.R, g.R', reshape(g.R(g.L==0), 1, []));
    
    % Compute coefficients as on RHS
    Iop.V = -Xi_coeff(K, R, Kp, Rp, K1, R1) .* double(g.L==g.L');
    
    % Locate i,j,z positions for non-zero elements
    ind = find(logical(Iop.V));
    [i, j, z] = ind2sub(size(K), ind);
    
    % Map to diagonal in energy, shifted to correct global index in I and J
    Iop.I = g.INeps + reshape(i-1, 1, [])*g.Neps;
    Iop.J = g.INeps + reshape(j-1, 1, [])*g.Neps;
    Iop.Z = repmat(reshape(zin.nubar(z), 1, []), [g.Neps, 1]);
    Iop.V = repmat(reshape(Iop.V(ind),   1, []), [g.Neps, 1]);
    
    if ~isempty(Iop)
        M = append_IJVZ(M, Iop, [], []);
    end

end

