function [fdot, b_z, bz_wtime] = ebe_solution(t, X, p, M, return_Jacobian, b_z_star)
    
    % New b_z
    [b_z, ~, bz_wtime] = new_bz(t, X, p, M);
    
    % Evaluate change in b_z for RHS, otherwise original b_z
    if ~isempty(b_z_star) && ~return_Jacobian
        b_z2 = b_z - b_z_star;
    else
        b_z2 = b_z;
    end
    
    % Create either dFdt or Jacobian based on flag
    if return_Jacobian
        fdot = sparse(M.m, M.n, M.Cpz*b_z2, M.grid.N, M.grid.N);
    else % fdot
        fdot = sparse(M.m, M.n, M.Cpz*b_z2, M.grid.N, M.grid.N) * X(:);
    end
    
end
