function [L1, L2] = L1_L2_coeff(L, Lp)
    % Note, R, RP is assumed to be carried in Delta coeff

    % Ensure working in double for arthimetic
    L = double(L);
    
    % L1 Factors
    L1_Lp1 = (L + 1.0) ./ (2.0 .* L + 3.0);
    L1_Lm1 =         L ./ (2.0 .* L - 1.0);
    
    % L2 Factors
    L2_Lp1 = (L + 2.0) ./ 2.0;
    L2_Lm1 = (1.0 - L) ./ 2.0;
    
    % Lp=L\pm1, using ints to be careful
    L = cast(L, 'int64');
    Lp = cast(Lp, 'int64');
    Lp1 = double(Lp == L+1);
    Lm1 = double(Lp == L-1);
    
    % Assemble L1 Term
    L1 = L1_Lp1 .* Lp1 + L1_Lm1 .* Lm1;
    L2 = L2_Lp1 .* Lp1 + L2_Lm1 .* Lm1;
    
end