
function [I, J, V, sid1, sid2] = convert_sparse_format(A, ind_bc)
    [I, J, V] = find(A); % Get the I/J/V coordinates as doubles
    [V, sid1] = sortrows([I, J, V], [1, 2]); % Sort by row, then column
    if nargin==2
        I = int32(V(:, 1)); % Store sorted J indices, matching V
        J = int32(V(:, 2)); % Store sorted J indices, matching V
        [ii, jj] = ind2sub(size(A), ind_bc);
        ii = int32(ii);
        jj = int32(jj);
        sid2 = zeros(numel(ii), 1);
        Nii = numel(ii);
        parfor k = 1:Nii
            [sid2(k)] = find( I == ii(k) & J == jj(k) );
            disp(k/Nii)
        end
    end
    I = int32(find(diff([0; int32(V(:, 1))]))); % compress row to pointers to first element
    J = int32(V(:, 2)); % Store sorted J indices, matching V
    V = V(:, 3);
    I(end+1) = nnz(V) + 1;
end
