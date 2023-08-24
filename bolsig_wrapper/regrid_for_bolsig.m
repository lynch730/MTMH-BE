

clear;clc;

NL_max = 950;

path_str = 'C:\Users\lynch\Documents\GitHub\boltzmann_solvers\lxcat_files\';
file_in_str = [path_str, 'Laporta_N2_vib_set.txt'];
file_out_str = [path_str, 'Laporta_N2_vib_set_new2.txt']; 
% 
% NL_max = 950;
% path_str = 'C:\Users\lynch\Desktop\MultiBolt-master (1)\MultiBolt-master\cross-sections\';
% file_in_str = [path_str, 'Laporta_N2_vib_set_new.txt'];
% file_out_str = [path_str, 'Laporta_N2_vib_set_new2.txt']; 

%%%%

lines = reshape(regexp(fileread(file_in_str), '\r\n|\r|\n', 'split'), [], 1);

% create empty fields
ind = reshape(find(startsWith(lines,  '------------')), 2, []);

ind_dash_all = true(size(lines));

for j = 1:size(ind, 2)

    ind_dash = (ind(1, j)+1):(ind(2, j)-1);

    data = lines(ind_dash, 1);
    data = reshape( sscanf(sprintf(' %s', data{:}), '%f %f \n', [1,Inf]), 2, [])';
    data = log10(data);
    Nr1 = size(data, 1);

    Nr2 =  min(NL_max, Nr1);
    x2 = linspace(data(1, 1), data(end, 1), Nr2);
    y2 = 10.0 .^ interp1(data(:, 1), data(:, 2), x2, 'pchip');
    y2(isnan(y2) | isinf(y2) | y2<0.0) = 0.0;

    lines( ind_dash(:), 1) = cell(Nr1, 1);
    lines( ind_dash(1:Nr2), 1) = strcat(sprintfc('%d',  10.0 .^ x2(:)), sprintfc('\t%d', y2(:)));

    ind_dash_all( ind_dash((Nr2+1):end) ) = false;
    
    fprintf('\n %i/%i', [j, size(ind, 2)]);
    
end

lines2 = lines(ind_dash_all);

fid = fopen(file_out_str,'w');
fprintf(fid,'%s\n',lines2{:});
fclose(fid);
