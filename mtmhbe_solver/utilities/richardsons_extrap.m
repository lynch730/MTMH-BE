function p = richardsons_extrap(x, y, nspace, yexact, ntrend)
    
    if ~exist('yexact', 'var')
        yexact = [];
    end
    if ~exist('ntrend', 'var')
        ntrend = [];
    end

    p.x = x(:);
    p.y = y(:);

    % Zero arrays
    N = numel(y);
    p.order_est = nan(N, 1);
    p.best_est = nan(N, 1);
    p.gci_fine = nan(N, 1);
    p.gci_course = nan(N, 1);
    p.asym_ratio = nan(N, 1);

    % For all cells greater than 3
    % i-2 is most coarse of 3
    % i is most fine
    assert(N>=3, 'Insufficient grid points for extrapolation')
    for i = N:-1:3
        order_est = log(abs((y(i-2)-y(i-1))/(y(i-1)-y(i))))/log(nspace);
        
        nexp = nspace.^order_est;
        best_est = y(i)+(y(i)-y(i-1))./(nexp-1.0);
        gci_fine =   1.25*abs((y(i)-y(i-1))/y(i))./(nexp-1.0);
        gci_course = 1.25*abs((y(i-1)-y(i-2))/y(i-1))./(nexp-1.0);
        asym_ratio = abs( gci_course ./ (nexp .* gci_fine) );
        if asym_ratio>1.1 || asym_ratio<0.9 || ~isreal(order_est)
            break
        else
            p.asym_ratio(i-2:i) = asym_ratio;
            p.gci_fine(i-2:i) = gci_fine;
            p.gci_course(i-2:i) = gci_course;
            p.best_est(i-2:i) = best_est;
            p.order_est(i-1:i) = order_est;
        end

    end

   if ~isempty(yexact)
        p.error_best = abs((y-yexact)./yexact);
   else
%         p.error_best = abs((y-p.best_est(end))./p.best_est(end));
        p.error_best = abs((y-p.y(end))./p.y(end));
   end

    p.in_domain_conv = abs(p.asym_ratio) < 1.1;
   

    if ~isempty(ntrend)
        x = x(ntrend:end);
        y = p.error_best(ntrend:end);
    else
        if ~any( p.in_domain_conv )
            return
        end
        x = x(p.in_domain_conv);
        y = p.error_best(p.in_domain_conv);
    end

    p.trend_line = polyfit(log10(x), log10(y), 1);

    p.xtrend = linspace(min(x), max(x), 100);
    p.ytrend = 10.0.^(polyval(p.trend_line, log10(p.xtrend)));
    p.order = p.trend_line(1);

end