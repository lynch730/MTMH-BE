%% Species Fractions
function [N_z, N0] = species_fractions(gas, zin, const)
    % X_total is the raw fraction of each species in the mixture
    % 

    % Constants
    if isempty(const)
        const = boltz_constants;
    end
    
    % Helpers
    X_u = reshape(gas.spec_frac, 1, []);
    Tref = gas.Tgas;
    if Tref<=0
        Tref = 300.0;
    end
    N0 = gas.press_Pa ./ (const.KB .* Tref);
    
    % Excitation temp 
    Texc_e = gas.Texc * const.KB / const.QE;
    
    % Start with Nk associated with each's multiplier's base species
    N_z = zeros(zin.Nz, 1);
    N_z(zin.coll) = N0 .* reshape(X_u( zin.u_z(zin.coll) ), [], 1);
    
    % Partition Function (size of total x-fractions)
    if zin.ISM == 0
        
        Z = zin.gratio.*exp(-zin.de./Texc_e);
        
        ind = logical(zin.de) & ~zin.is_upper_state & ~zin.is_thermal;
        N_z(ind) = N_z(ind) ./ (1.0 + Z(ind));
        
        ind = logical(zin.de) & zin.is_upper_state & ~zin.is_thermal;
        N_z(ind) = N_z(ind) .* Z(ind) ./ (1.0 + Z(ind));
        
    elseif zin.ISM == 1 % Standard approach for Ensemble
        
        % Identify Upper State, so that ground states get correct zero de
        % Note, de stored value is not zero for ground state, because
        %    we use it to ID the Z index as inelastic (below)
        isup = zin.is_upper_state;
        is_scaled = logical(zin.de); % only apply to inelastic terms

        % Partition functions
        % Technically, we should ensure gratio of root z species is always set
        % to one, as it's value is carried in gratio, but likely rare to occur
        % (Note: zero gratio prevents this from apply to total species!)
        if Texc_e>0
            Z_z = zin.gratio.*exp(-zin.de .* double(isup & is_scaled) ./ Texc_e);
        else
            Z_z = zeros(zin.Nz, 1);
            Z_z(~isup & is_scaled) = 1.0;
        end

        % Sum Partition functions associated with the same root species
        %   (+1 so no zero index accumulation from terms without u_z)
        Z_u = accumarray(zin.u_z+1, Z_z); % sum on 1+u species index
        Z_u = Z_u(2:end); % Don't use the +1 garbage bin

        % For all inelastic/superelastic, divide by total species Z
        % and multiply by ground state
        N_z(is_scaled) = N_z(is_scaled) .* Z_z(is_scaled) ./ Z_u( zin.u_z(is_scaled) ) ;
        
    elseif zin.ISM == 2 % Standard approach for Ensemble
        
    end
    
end
