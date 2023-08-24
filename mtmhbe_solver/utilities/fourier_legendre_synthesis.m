function [v, phase_base] = fourier_legendre_synthesis( p, M, vec, new_time, is_mom)

    if ~exist('is_mom', 'var')
        is_mom = 0;
    end

    if ~exist('new_time', 'var') 
        new_time = [];
    end

    g = M.grid;
    if is_mom
        L0 = find(g.L==1);
        NL0 = numel(L0);
    else
        L0 = find(g.L==0);
        NL0 = numel(L0);
    end

    if isempty(new_time)

        %% Steady state solution (final value)
        K = reshape(double(g.K(L0)), [1, 1, NL0]);
        R = reshape(double(g.R(L0)), [1, 1, NL0]);
        
        phase_base = linspace(0, 2.0*pi, 200);
        cs = cos(phase_base .* K) .* (1.0-R) + sin(phase_base .* K) .* R;
    
        vec2 = vec;
        if min(size(vec)) > 1
            vec2 = reshape(vec2(:, :), [1, size(vec2)]);
        else
            vec2 = reshape(vec2, 1, 1, []);
            assert(NL0==size(vec2,3), 'Vector input must have dimension NL0')
        end
        
        v = zeros(size(vec2, 1), numel(phase_base));
    
        for j = 1:size(vec2, 1)
            for ik = 1:NL0
                v(j, :) = v(j, :) +  vec2(j, 1, ik) .* cs(1, :, ik);
            end
        end

    else
        
        %% Map to user-given


        sz = size(vec);
        if ndims(vec) == 3

            K = reshape(double(g.K(L0)), [1, NL0]);
            R = reshape(double(g.R(L0)), [1, NL0]);

            new_time = reshape(new_time, 1, 1, []);
            phase_base = p.omega .* new_time;
            cs = cos(phase_base .* K) .* (1.0-R) + sin(phase_base .* K) .* R;
    
            v = zeros(sz(1), numel(new_time));
            for j = 1:sz(1)
                for ik = 1:NL0
                    mean_loc = interp1(p.time.array, ...
                                       reshape(vec(j, ik, :), 1, []), ...
                                       new_time, 'linear');
%                     mean_loc = reshape(mean_loc, 1, []);
                    v(j, :) = v(j, :) + reshape(mean_loc.*cs(1, ik, :), 1, []);
                end
            end

        else


            K = reshape(double(g.K(L0)), [NL0, 1]);
            R = reshape(double(g.R(L0)), [NL0, 1]);

            new_time = reshape(new_time, 1, []);
            phase_base = p.omega .* new_time;
            cs = cos(phase_base .* K) .* (1.0-R) + sin(phase_base .* K) .* R;
            
            v = zeros(1, numel(new_time));
            for ik = 1:NL0
                mean_loc = interp1(p.time.array, ...
                                   reshape(vec(ik, :), 1, []), ...
                                   new_time, 'linear');
                v(1, :) = v(1, :) + mean_loc .* cs(ik, :);
            end            
        end


    end

end