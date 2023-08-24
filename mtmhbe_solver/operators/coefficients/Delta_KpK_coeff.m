function dkpk = Delta_KpK_coeff(K, R, Kp, Rp)
    
    % Ensure working in integers
    K = cast(K, 'int64');
    R = cast(R, 'int64');
    Kp = cast(Kp, 'int64');
    Rp = cast(Rp, 'int64');
    
    % Delta Kp, K factors 
    dkA = double(Kp==(K-1));
    dkB = double(Kp==(K+1));
    dkC = double(Kp==0) .* double(K==1);
    
    % Rmatch
    rmatch = double(Rp==R);
    
    % See Scale Factor
    dkpk = 0.5.*(dkA + dkB + dkC).*rmatch;
    
end