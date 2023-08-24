function grid = energy_grid(grid)
    
    %% Step 1) LKR sub-matrices
        
        % Set the LKR sub-matrices
        [L, K, R] = lkr_indices(grid.NL, grid.NK, grid.FL_order);
        
        % Set size of LKR set and major dimension N
        grid.NLK = numel(L);
        
        % Store LKR
        grid.L = L;
        grid.K = K;
        grid.R = R;

        % Custom Field
        if strcmp(grid.grid_case, 'custom')
            assert(isfield(grid, 'EC'), 'EC must be defined')
            grid.EC = unique(grid.EC);
            grid.Neps = numel(grid.EC);
        end

        % Store Temp
        Neps = grid.Neps;

        % Grid index arrays
        grid.INeps = reshape(1:Neps, [], 1);
        
        % Selectors for isotropic terms
        grid.L_is_0 = find(L==0);
        grid.NL0 = numel(grid.L_is_0);
        grid.NL1 = numel(find(grid.L==1));
        grid.ind_L0 = grid.INeps + reshape( (grid.L_is_0-1)*grid.Neps, 1, []);
        
        % Sizes for each type
        grid.N = grid.NLK * Neps;
        
        % Even Terms
        grid.L_is_even = ~mod(L, 2);
        
    
    %% Step N) Tiling matrices
    
        % % L=0 Tile
        % im = grid.L_is_0;
        % grid.tile.L0 = reshape((im-1)*Neps, 1, 1, []);
        % 
        % % Even L>0 Tile
        % im = find(grid.L>0 & grid.L_is_even);
        % grid.tile.Le = reshape((im-1)*Neps, 1, 1, []);
        % 
        % % Odd L>0 Tile
        % im = find(grid.L>0 & ~grid.L_is_even);
        % grid.tile.Lo = reshape((im-1)*Neps, 1, 1, []);
        
    %% Step 2) Define Odd/Even Energy Grids
        
        % Define grid at cell walls
        if strcmp(grid.grid_case, 'custom')  
            grid.EC = reshape(grid.EC, [], 1);
            grid.ER = grid.EC + 0.5.*[diff(grid.EC); grid.EC(end)-grid.EC(end-1)];
            grid.eV_min = grid.ER(1);
            grid.eV_max = grid.ER(end);
            grid.EC = [];
        else
            grid.ER = grid_spacing(grid.eV_min, grid.eV_max, ...
                                     Neps, grid.grid_case);
        end
        
        % Even grid
        grid.EL = [0; grid.ER(1:end-1)];
        grid.Evol_0 = (2.0./3.0) .* (grid.ER.^1.5 - grid.EL.^1.5); % e^1/2
        grid.Evol_1 = (1.0./2.0) .* (grid.ER.^2.0 - grid.EL.^2.0); % 
        grid.Evol_2 = (2.0./5.0) .* (grid.ER.^2.5 - grid.EL.^2.5);
        grid.Evol_3 = (1.0./3.0) .* (grid.ER.^3.0 - grid.EL.^3.0);
        grid.Evol_4 = (2.0./7.0) .* (grid.ER.^3.5 - grid.EL.^3.5);
        grid.EC = grid.Evol_2./grid.Evol_0;

        assert(~any(diff(grid.EC)<0), 'Grid not sorted');
        assert(~any(isnan(grid.ER)), 'NaNs in custom array!');
        assert(~any(isnan(grid.EC)), 'NaNs in custom array!');

        % Odd grid
        grid.OL = 0.5.*(grid.ER + grid.EL);
        grid.OL(1) = 0.0;
        grid.OR = [grid.OL(2:end); 2*grid.OL(end)-grid.OL(end-1)];
        grid.Ovol_0 = (2.0./3.0) .* (grid.OR.^1.5 - grid.OL.^1.5);
        grid.Ovol_1 = (1.0./2.0) .* (grid.OR.^2.0 - grid.OL.^2.0);
        grid.Ovol_2 = (2.0./5.0) .* (grid.OR.^2.5 - grid.OL.^2.5);
        grid.Ovol_3 = (1.0./3.0) .* (grid.OR.^3.0 - grid.OL.^3.0);
        grid.Ovol_4 = (2.0./7.0) .* (grid.OR.^3.5 - grid.OL.^3.5);
        grid.OC = grid.Ovol_2./grid.Ovol_0;
                
        % Finite volume of each cell, 
        %    _0 over mass
        %    _1 over momentum
        %    _2 over energy
        % grid.EO_vol_0 = reshape(grid.Evol_0.*double(grid.L_is_even) ...
        %                 + grid.Ovol_0.*double(1.0-grid.L_is_even), [], 1);
        % grid.EO_vol_1 = reshape(grid.Evol_1.*double(grid.L_is_even) ...
        %                 + grid.Ovol_1.*double(1.0-grid.L_is_even), [], 1);
        % grid.EO_vol_2 = reshape(grid.Evol_2.*double(grid.L_is_even) ...
        %                 + grid.Ovol_2.*double(1.0-grid.L_is_even), [], 1);      
         
        % Left and right cell energies
        % Note that 1st even grid uses forward diff (1, 2) to estimate slope
        %           last even grid uses extrapolated OR(end), assuming zero
        %           value
        %           1st Ogrid left uses zero value at e=0
        %           last Ogrid right uses backward slope
        grid.ELC = [grid.EC(1); grid.EC(1:end-1)];
        grid.ERC = [grid.EC(2:end); grid.OR(end)];
        grid.OLC = [grid.EL(1); grid.OC(1:end-1)];
        OCNP = grid.OC(end) + 0.5.*diff(grid.OC(end-1:end));
        grid.ORC = [grid.OC(2:end); OCNP];

        % Mass Integration
        grid.E_INT_MASS = integral_matrix(grid, grid.Evol_0, grid.Evol_2, true);
        grid.O_INT_MASS = integral_matrix(grid, grid.Ovol_0, grid.Ovol_2, false);
        % grid.EO_INT_MASS = reshape(grid.E_INT_MASS.*double(grid.L_is_even) ...
        %                          + grid.O_INT_MASS.*double(1.0-grid.L_is_even), [], 1);
        grid.E_INT_MASS = reshape(grid.E_INT_MASS, 1, []);
        grid.O_INT_MASS = reshape(grid.O_INT_MASS, 1, []);

        % Momentum Integration
        grid.E_INT_MOM = integral_matrix(grid, grid.Evol_1, grid.Evol_3, true);
        grid.O_INT_MOM = integral_matrix(grid, grid.Ovol_1, grid.Ovol_3, false);
        % grid.EO_INT_MOM = reshape(grid.E_INT_MOM.*double(grid.L_is_even) ...
        %                          + grid.O_INT_MOM.*double(1.0-grid.L_is_even), [], 1);
        grid.E_INT_MOM = reshape(grid.E_INT_MOM, 1, []);
        grid.O_INT_MOM = reshape(grid.O_INT_MOM, 1, []);
        
        % Integration Kernel
        grid.E_INT_ENERGY = integral_matrix(grid, grid.Evol_2, grid.Evol_4, true);
        grid.O_INT_ENERGY = integral_matrix(grid, grid.Ovol_2, grid.Ovol_4, false);
        % grid.EO_INT_ENERGY = reshape(grid.E_INT_ENERGY.*double(grid.L_is_even) ...
        %                          + grid.O_INT_ENERGY.*double(1.0-grid.L_is_even), [], 1);
        grid.E_INT_ENERGY = reshape(grid.E_INT_ENERGY, 1, []);
        grid.O_INT_ENERGY = reshape(grid.O_INT_ENERGY, 1, []);

        % Grid helpers (for collision terms, mainly)
        [grid.ue, grid.uo] = grid_helpers(grid);

end


function y = integral_matrix(g, I1, I2, is_even)
    if is_even
        URC = g.ERC;
        UC = g.EC;
        ULC = g.ELC;
    else
        URC = g.ORC;
        UC = g.OC;
        ULC = g.OLC;
    end
    de = URC - ULC;
    de(1) = URC(1) - UC(1);
    de(end) = UC(end) - ULC(end);
    I3 = (I2 - UC .* I1)./de;
    V = [-I3(:), I1(:), I3(:)];
    J = [g.INeps-1, g.INeps, g.INeps+1];
    J(1, 1) = 1;
    J(end, 3) = g.Neps;
    y = accumarray(J(:), V(:), [g.Neps, 1]);
end

%% Grid spacing
function x = grid_spacing( emin, emax, N, grid_case, k)
    assert(emin>0.0, 'Minimum energy must be greater than zero!');
    assert(emin<emax, 'Minimum must be less than maximum')
    assert(N>=6, 'Grid should not be smaller than 10 cells')
    switch grid_case
        case 'log'
            x = reshape(10.^linspace(log10(emin), log10(emax), N), [], 1);
        case 'quad'
            if ~exist('k', 'var')
                k = 2;
            end
            x = (linspace(emin.^(1.0./k), emax.^(1.0./k), N)).^k;
        case 'linear'
            x = linspace(emin, emax, N);
        otherwise
            error('Bad value for grid_case, must be: log, quad, linear')
    end
    x = reshape(x, [], 1);
end


% Creates helpers grids
% In future, build these naturall!
function [ue, uo] = grid_helpers(g)
    
    % Even grid
    ue.Neps = g.Neps;
    ue.INeps = g.INeps;
    ue.is_even = true;
    ue.cell_volume = g.Evol_0;
    ue.eps_c_i = g.EC;
    ue.eps_l_i = g.EL;
    ue.eps_r_i = g.ER;
    ue.eps_i = [ue.eps_l_i; ue.eps_r_i(end)];
    ue.eps_c_i_diff2 = g.ERC - g.ELC;
    ue.bc_end_weight = 1.0 - ue.eps_c_i_diff2(end)./...
                            (ue.eps_r_i(end) - g.ELC(end));
    ue.ediff = (g.EC - g.ELC); % Nth left cell edge
    ue.ediff(1) = 1; % Arbitrary non-zero i=1 value
    ue.ediff = [ue.ediff; 1]; % Arbitrary non-zero N+1 value 

    % Odd grid
    uo.Neps = g.Neps;
    uo.INeps = g.INeps;
    uo.cell_volume = g.Ovol_0;
    uo.is_even = false;
    uo.eps_c_i = g.OC;
    uo.eps_l_i = g.OL;
    uo.eps_r_i = g.OR;
    uo.eps_i = [uo.eps_l_i; uo.eps_r_i(end)];
    uo.eps_c_i_diff2 = g.ORC - g.OLC;
    uo.bc_end_weight = 1.0 - uo.eps_c_i_diff2(end)./...
                            (uo.eps_r_i(end) - g.OLC(end));
    uo.ediff = (g.OC - g.OLC); % Nth left cell edge
    uo.ediff(1) = 1; % Arbitrary non-zero i=1 value
    uo.ediff = [uo.ediff; 1]; % Arbitrary non-zero N+1 value 

end
