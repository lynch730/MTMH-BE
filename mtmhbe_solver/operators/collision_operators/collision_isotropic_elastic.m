function [C_el_I, C_el_T] = collision_isotropic_elastic(proc, g, zin_thermal)
    
    %% Cell Eedge Values
    const = boltz_constants;

    % Right and left function value constructed from linear interpolation
    % upwind. Left cell edge takes eps_i and eps_i+1 (with weight wi)
    % Right cell edge takes eps_i+1 and eps_i+2

    % Right Edge Cell: eps^2 * sigma(eps)
    %  Right Cell Edge Flux = 0;
    Ai = extended_linear_interp(proc, g.ER);
    Ai = Ai .* g.ER .* g.ER;
    Aiend = Ai(end);
    Ai(end) = 0; % Right Edge BC
    
    % Left Cell Edge
    %  Boundary Condition: Left Edge Flux = 0
    Aim = [0; Ai(1:end-1)];
    
    %% Weights for Elastic Cooling
    
    % Weights
    wi = (g.ERC - g.EL) ./ (g.ERC - g.EC);
    
    % Left face of last cell is calculated from f(eps_N) and condition
    % f(eps_r,N)=0. Interpolation back to left face comes out as scaling
    % f(eps_N) by some factor
    wi(end) = (g.ER(end) - g.EL(end))./(g.ER(end)-g.EC(end));
    
    % Right cell with zero (arbitrary due to Ai BC above), but zero anyway
    wip = [wi(2:end); 0];
    
    %% Non-thermal Elastic Cooling

    % Assemble Coefficients
    C_el_I.V = [-Aim .* wi, ...
                 Ai .* wip - Aim .* (1.0 - wi), ...
                 Ai .* (1.0 - wip)];
    
    % TEMPORARY
    C_el_I.V(end, 2) = 0.0;
    C_el_I.V(end-1, 3) = 0.0;
    
    % Divide by cell volume
    C_el_I.V = C_el_I.V .* (2.0.*proc.mratio) .* const.GAMMA;
    C_el_I.V = C_el_I.V ./ g.Evol_0;
    
    % I/J Pointers
    C_el_I.I = repmat(g.INeps, [1, 3]);
    C_el_I.J = [g.INeps, g.INeps+1, g.INeps+2];
    C_el_I.Z = C_el_I.I*0 + proc.z_s;
    
    % Adjust J Pointers to be in bounds (arbitrary, associated values are zero)
    C_el_I.J(end-1, 3) = C_el_I.J(end-1, 2); % Sets constant V for 2nd to last cell
    C_el_I.J(end, 2:3) = C_el_I.J(end,   1); % Sets constant V for last cell
    
    
    %% Thermal terms (kbT included in beta)
    
    % Recycle Elastic Terms, Divide by (eps_i+1 - eps_i)
    Ai = Ai ./ (g.ERC -  g.EC);
    Aim = [0; Ai(1:end-1)];
    
    % Collect Terms
    C_el_T.V = [Aim, -Ai - Aim, Ai];
    C_el_T.V(end, 2) = C_el_T.V(end, 2) - Aiend ./ (g.ER(end) - g.EC(end));

    % I/J Pointers
    C_el_T.I = repmat(g.INeps, [1, 3]);
    C_el_T.J = [g.INeps-1, g.INeps, g.INeps+1 ];
    C_el_T.Z = C_el_T.I*0 + zin_thermal;
    
    % Adjust J Pointers to be in bounds (arbitrary, associated values are zero)
    C_el_T.J(1, 1) = 1; 
    C_el_T.J(end, 3) = g.Neps; 
    
    % Scale and Divide by cell volume
    C_el_T.V = C_el_T.V .* (2.0.*proc.mratio) .* const.GAMMA;
    C_el_T.V = C_el_T.V ./ g.Evol_0;
    
end