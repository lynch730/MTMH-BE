function sig = extended_linear_interp(proc, xq)
    
    % Apply linear approx. 
    if proc.log_interp_flag
        xq(xq<=0) = 1e-20;
        sig = 10.0.^proc.sigma_log(log10(xq));
        sig(isnan(sig)) = 0.0;
    else
        sig = proc.sigma_linear(xq);
    end

    % If extrapolation is on and process does not exented to requested
    % grid, apply extrapolation
    is_oob = xq > proc.eps_j(end);
    if any(is_oob)
        sig(is_oob) = proc.sigma_j(end);
        if proc.extrap
            sig(is_oob) = sig(is_oob) .* log10(xq(is_oob)) ./ xq(is_oob);
            sig(is_oob) = sig(is_oob) ./ (log10(proc.eps_j(end))./ proc.eps_j(end));
        end
    end
    
    % Ensure no nans
    sig(isnan(sig)) = 0.0;
    
end
