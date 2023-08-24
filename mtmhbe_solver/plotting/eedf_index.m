
function [iplot, Np] = eedf_index(M, pcase)
    
    % Set
    % Helpers
    g = M.grid;

    % Select distros to plot
    switch pcase
        case 'isotropic'
            iplot = find(g.L==0);
        case 'isotropic-real'
            iplot = find(g.L==0 & g.R==0);
        case 'isotropic-dc'
            iplot = find(g.L==0 & g.K==0);
        otherwise
            iplot = 1:g.NLK;
    end

    Np = numel(iplot); % number of distros to plot

end
