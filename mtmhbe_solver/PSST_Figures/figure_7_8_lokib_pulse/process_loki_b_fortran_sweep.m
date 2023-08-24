

clear; clc;

fname = 'fort_Ne_sweep_ebar';
str = fileread([fname, '.txt']);
lines = regexp(str, '\r\n|\r|\n', 'split');
lines(cellfun('isempty',lines)) = [];
ind = find(strcmp(lines, ' ------------------------------------------' ));

Ndat = numel(lines(1:ind(1)-2));
dat = zeros(Ndat, 3);
for i = 1:Ndat
    line_loc = split(strtrim(lines{i}));
    dat(i,:) = str2double(line_loc);
end

Ntimes = ind(2) - ind(1) -1;
wall_time = zeros(Ntimes, 1);
for i = 1:Ntimes
    line_loc = strsplit(lines{ind(1)+i}, ',');
    wall_time(i) = str2double(line_loc{2});
end

times = reshape(dat(:, 1), [], Ntimes);
times = times(:, 1);

Nt = numel(times);
wall_time = wall_time * Nt / 1000.0;

ebar = reshape(dat(:, 3), [], Ntimes);

fname = 'time_dependent_fortran_Neps_sweep.mat';
save(fullfile('PSST_Figures','figure_7_8_lokib_pulse',fname), 'times', 'ebar', 'wall_time')

