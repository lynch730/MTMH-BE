function mset = append_IJVZ(mset, oper, i_tile, j_tile)
    
    % Convert integer types and flatten
    I = cast(oper.I(:), mset.int_type);
    J = cast(oper.J(:), mset.int_type);
    Z = cast(oper.Z(:), mset.int_type);
    V = oper.V(:);
    
    % Clear any zero entries in V
    ind = find(~logical(V));
    if ~isempty(ind)
        I(ind) = [];
        J(ind) = [];
        V(ind) = [];
        Z(ind) = [];
    end
    
    % Apply tiling if given
    if ~isempty(i_tile) && ~isempty(j_tile)

        % Tiles as columns
        i_tile = reshape(i_tile, 1, []);
        j_tile = reshape(j_tile, 1, []);

        % Row-Column scaling of I and J indices, reshaped as single column
        I = reshape(I + cast((i_tile-1)*mset.Neps, mset.int_type), [], 1);
        J = reshape(J + cast((j_tile-1)*mset.Neps, mset.int_type), [], 1);

        % Equivalent copying of Z and V, with no change
        Z = reshape(repmat(Z, 1, numel(i_tile)), [], 1);
        V = reshape(repmat(V, 1, numel(i_tile)), [], 1);

    end

    % Check Rows for OOB Indices
    ind_oob = find( I(:) > mset.N |  I(:) < 1 | J(:) > mset.N | J(:) < 1 );
    if ~isempty(ind_oob)
        fprintf('\nOut of bounds found in I/J, coordinates:')
        for k = 1:numel(ind_oob)
            fprintf('\n\t %i %i %i', ...
                    [I(ind_oob(k)), J(ind_oob(k)), Z(ind_oob(k))] );
        end
        error('I/J OOB Error')
    end
    
    % Append to stack 
    %   (If this takes a long time, the cause is re-sizing
    %    fields. Solution is to adjust Nnz_est in Matrix main)
    k = mset.row_cnt;
    if ~isempty(I)
        k = k:(k+numel(I)-1);
        mset.I(k) = I(:);
        mset.J(k) = J(:);
        mset.Z(k) = Z(:);
        mset.V(k) = V(:);
        mset.row_cnt = k(end)+1;
    end

end

