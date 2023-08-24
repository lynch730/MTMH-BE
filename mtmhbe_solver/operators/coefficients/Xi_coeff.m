function Xi = Xi_coeff(K, R, Kp, Rp, K1, R1)
    
    % Ensure working in integers
    % We have to be very very careful to use correct types!!!!
    K = cast(K, 'int64');
    R = cast(R, 'int64');
    Kp = cast(Kp, 'int64');
    Rp = cast(Rp, 'int64');
    K1 = cast(K1, 'int64');
    R1 = cast(R1, 'int64');
    
    % K Kronecker Deltas
    KdelA = double(K == (K1 + Kp));
    KdelB = double(K == (K1 - Kp));
    KdelC = double(K == (Kp - K1));
    
    % R Kronecker Deltas
    RdelA = 1.0 - 2.0 .* double(R==0) .* double(R1==1);
    RdelB = 1.0 - 2.0 .* double(R==1) .* double(R1==0);
    RdelC = 1.0 - 2.0 .* double(R==1) .* double(R1==1);
    
    % Sum K and R's
    Xi = KdelA .* RdelA + KdelB .* RdelB + KdelC .* RdelC;

    % Add Coupled from Neps_00 to all DC terms
    Xi = Xi + double(K==0).*double(K1==0).*double(Kp==0);

    % Coupling on R's
    Xi = Xi .* double(Rp==round(abs(R1-R)));

    % Divide by 2, extra 2 for K=0
    Xi = Xi ./ (2.0 .*(1.0 + double(K==0)));
    
end