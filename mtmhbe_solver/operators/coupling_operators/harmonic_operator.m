function M = harmonic_operator(M, g, zin)
    
    % Location same L terms, Same K, Opposite real/imaginary
    [i, j] = find(g.L'==g.L & g.K'==g.K & g.R'~=g.R); 
    
    % I/J/Z main diagonals
    H.I = g.INeps + reshape(g.Neps*(i-1), 1, []);
    H.J = g.INeps + reshape(g.Neps*(j-1), 1, []);
    H.Z = H.I*0 + zin.omega;
    
    % Value from K and R as doubles
    H.V = -double(g.K(i)) .* (-1).^double(g.R(i)); % -K*(-1)^R
    
    % Reshape each sub-matrix instance to columns and copy down Neps (diag)
    H.V = repmat(reshape(H.V, 1, []), g.Neps, 1);
    
    % Append to matrix stack, no tiling
    M = append_IJVZ(M, H, [], []);
    
end


    

