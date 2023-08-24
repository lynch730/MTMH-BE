
clear, clc;

bd = pulse_settings;

base_name = 'loki_b_pulse';
dname = fullfile('data','loki_b_pulse');

bd.grid.Neps = 200;
NK_array = 5:5:50;
NL_array = 5:5:50;

for i = 1:numel(NK_array)   
%     for j = 1:numel(NL_array)
        
        bd.grid.NK = NK_array(i);
        bd.grid.NL = NL_array(i);

        M = matrix_main(bd.xsec, bd.grid, bd.paths);
        
        fname = [base_name, '_N', sprintf('%i', bd.grid.Neps), ...
                            '_NL', sprintf('%i', bd.grid.NL), ...
                            '_NK', sprintf('%i', bd.grid.NK), ...
                            '.mat']; % For grid study
        
        export_operators(M, bd.paths, fullfile(dname, fname));
        
        fprintf('\n %i', i)
    
%     end
end
