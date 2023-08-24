function [M, xsec, zin] = collision_operators(M, g, zin, xsec)

    % Diagonal tiling arrays
    isotropic = find( g.L(:)==0 );
    aniso_even = find( g.L(:)>0 & ~mod(g.L(:), 2) );
    aniso_odd  = find( g.L(:)>0 &  mod(g.L(:), 2) );
    
    % Loop Processes
    for i = 1:xsec.Ns
        
        % Compute next z-index for process
        [zin, xsec.proc(i), zTherm] = z_index_collision(zin, xsec.proc(i), xsec.ISM);
        
        % Compute collision operators on odd and even grids
        %  Returns source for even to use on isotropic grids
        [sink.odd, ~, u.odd] = collision_process(xsec.proc(i), g.uo);
        [sink.even, source, u.even] = collision_process(xsec.proc(i), g.ue);
        
        % Anisotropies
        if ~isempty(aniso_even)
            M = append_IJVZ(M, sink.even, aniso_even, aniso_even);
        end
        if ~isempty(aniso_odd)
            M = append_IJVZ(M, sink.odd, aniso_odd, aniso_odd);
        end
        
        % For Elastic Collisions Isotropic Terms
        if xsec.proc(i).is_elastic 
            [C_el_I, C_el_T] = collision_isotropic_elastic(xsec.proc(i), g, zTherm);
            M = append_IJVZ(M, C_el_I, isotropic, isotropic);
            M = append_IJVZ(M, C_el_T, isotropic, isotropic);
        else % Normal Inelastic, applied to isotropic
            M = append_IJVZ(M, sink.even, isotropic, isotropic);
            if ~isempty(source)
                M = append_IJVZ(M, source, isotropic, isotropic);
            end
        end
        
        % Store Rate Integral for Even Terms 
        % Make Negative, to cancel the negative given as a sink
        M.rates.integral(i, :) = -sink.even.integral;
        M.rates.zid(i) = xsec.proc(i).z_s;

    end
    
    % Non-Conservative rates for nubar eval
    ind_nonc =  [xsec.proc(:).is_nonconservative]';
    zin.nubar_eval.any_nonc = any(ind_nonc);
    if  zin.nubar_eval.any_nonc
        zin.nubar_eval.z_s = [xsec.proc(ind_nonc).z_s]';
        zin.nubar_eval.ind_att = [xsec.proc(ind_nonc).is_attachment];
        zin.nubar_eval.integral = M.rates.integral(ind_nonc, :);
    end
    
end




    