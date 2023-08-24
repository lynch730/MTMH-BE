function [zin, proc, zin_thermal] = z_index_collision(zin, proc, ISM)
        
    % Default Thermal Index to Empty
    zin_thermal = [];

    % Switch by ISM Type
    if ISM == 2
        
        % Boolean Tests for paired ground state species Note that ISM-2
        % cross_section_main is different than for ISM-0/1 as it uses the
        % product species and u index directly, rather than as a
        % placeholder to identify the root species. Therefore if logic in
        % that subroutine is correct, every process s will have an
        % associated total species here, and tripping assert below
        % indicates a real problem with the code
        zcoll = find(zin.u_z == proc.u_s);
        assert(numel(zcoll)==1, 'Multiple valid z indices found!');
        [zin, proc] = add_process_to_existing(zin, proc, zcoll);
        
        % Add Thermal-Elastic Term
        if proc.is_elastic || proc.is_effective
            [zin, proc] = increment_z_coll(zin, proc, true);
            zin_thermal = zin.Nz;
        end
    
    else
        
        if proc.is_excitation
            
            % Super-elastics always increment, as do ISM-0
            if proc.is_superelastic || ISM==0
                [zin, proc] = increment_z_coll(zin, proc, false);
            else % Try to find ground state to use
                
                % Boolean Tests for paired ground state species 
                % looking for a previous inelastic (not superelastic)
                % with same root species, use de to preclude the total
                % states and other z parameters
                z_with_same_u = zin.u_z == proc.u_s;
                z_is_excit = zin.de ~= 0;
                z_is_ground = ~zin.is_upper_state;
                
                % Search for any z with an
                zcoll = find(z_with_same_u & z_is_excit & z_is_ground);
                if ~isempty(zcoll)
                    assert(numel(zcoll)==1, 'Multiple valid z indices found!');
                    [zin, proc] = add_process_to_existing(zin, proc, zcoll);
                else % Increment z
                    [zin, proc] = increment_z_coll(zin, proc, false);
                end

            end

        else % Elastic, effective, ionization, or attachment, using total species
            
            % Non-Thermal Piece from total species
            zcoll = (zin.coll(1)-1) + proc.u_s; % Locate in user species  
            [zin, proc] = add_process_to_existing(zin, proc, zcoll);
            
            % If Elastic, add thermal component
            if proc.is_elastic || proc.is_effective
                [zin, proc] = increment_z_coll(zin, proc, true);
                zin_thermal = zin.Nz;
            end
            
        end
        
    end
    
end

% Add to existing z index
function [zin, proc] = add_process_to_existing(zin, proc, zcoll)
    proc.z_s = zcoll;
    zin.s_z{proc.z_s} = [zin.s_z{proc.z_s}; proc.s]; 
end

%% Create new Z index
function [zin, proc] = increment_z_coll(zin, proc, is_thermal)
    
    % New ground state, increment
    zin.Nz = zin.Nz + 1;
    
    % Fill Index Arrays
    zin.coll(end+1) = zin.Nz;
    zin.u_z(end+1) = proc.u_s;
    zin.is_upper_state(end+1) = proc.is_superelastic;
    zin.is_thermal(end+1) = is_thermal;
    zin.s_z{end+1} = proc.s;
    
    % Transition energy, set to zero for thermal
    if is_thermal
        zin.de(end+1) = 0;
        zin.gratio(end+1) = 0;
    else
        zin.de(end+1) = proc.de;
        zin.gratio(end+1) = proc.gratio;
        proc.z_s = zin.Nz;  % Save z for this process, but not thermal
    end
    
end