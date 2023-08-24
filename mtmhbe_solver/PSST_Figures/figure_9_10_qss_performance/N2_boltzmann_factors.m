


%% BOltzmann statistics, assuming all g=1
function X = N2_boltzmann_factors(de, Texc)
    % Convert Texc to eV
    c = boltz_constants;
    Texc = Texc * c.KB ./ c.QE;
    if abs(Texc) > 100*eps
        Z = exp(-de./Texc);
        Zsum = sum(Z(:));
        X = Z./Zsum;
    else
        X = zeros(size(de));
        X(1) = 1.0;
    end
end