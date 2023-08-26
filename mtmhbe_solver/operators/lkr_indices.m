function [L, K, R] = lkr_indices(NL, NK, term_order)
 
    % Base Arrays
    L = 0:NL-1;
    K = 0:NK-1;
    R = [0, 1];
    
    % Choose ordering, fastest to slowest change left to right
    switch term_order 
        case 1
            [R, K, L] = ndgrid_fast(R, K, L); % Probably works?
        case 2 
            [K, R, L] = ndgrid_fast(K, R, L); % Works
        otherwise
            error('Bad FL Input')
    end
        
    % Extract the valid components of L/K to return
    if NK>1 % AC Field
        ind = mod(L + K, 2) == 0 & ~(K==0 & R==1);
    else % DC Field
        ind = R==0;
    end
    
    L = L(ind);
    K = K(ind);
    R = R(ind);

end

%% Faster form of ndgrid
function [A, B, C] = ndgrid_fast(a, b, c)
    s = [numel(a), numel(b), numel(c)];
    A = repmat(reshape(a, [s(1), 1, 1]), [1, s(2), s(3)]);
    B = repmat(reshape(b, [1, s(2), 1]), [s(1), 1, s(3)]);
    C = repmat(reshape(c, [1, 1, s(3)]), [s(1), s(2), 1]);
end
