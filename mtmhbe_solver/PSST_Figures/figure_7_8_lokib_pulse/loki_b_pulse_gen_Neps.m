
clear, clc;

bd = pulse_settings;

base_name = 'loki_b_pulse';
dname = fullfile('data','loki_b_pulse');

bd.grid.NK = 8;
bd.grid.NL = 8;

Neps_array = [50:50:500, 5000];

for i = 1:numel(Neps_array)   

    bd.grid.Neps = Neps_array(i);

    M = matrix_main(bd.xsec, bd.grid, bd.paths);
    
    fname = [base_name, '_N', sprintf('%i', bd.grid.Neps), ...
                        '_NL', sprintf('%i', bd.grid.NL), ...
                        '_NK', sprintf('%i', bd.grid.NK), ...
                        '.mat']; % For grid study
    
    export_operators(M, bd.paths, fullfile(dname, fname));
    
    fprintf('\n %i', i)

end
