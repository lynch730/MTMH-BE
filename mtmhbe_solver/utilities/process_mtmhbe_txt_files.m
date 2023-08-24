
function data = process_mtmhbe_txt_files(varargin)
    assert(nargin>=1, 'must have at least one input!');
    data = convert_data_file(varargin{1});
    for i = 2:nargin
        new_data = convert_data_file(varargin{i});
        data = pair_data(data, new_data);
    end
end

function data = convert_data_file(fname)

    lines = regexp(fileread(fname), '\r\n|\r|\n', 'split')';
    lines(cellfun('isempty',lines)) = [];
    ind = find(strcmp(lines, ' ------------------------------------------' ));
    
    ind_ebar_start = find(strcmp(lines, ' Mean Energies' )) + 1;
    ind2 = ind_ebar_start:(ind-2);
    ebar_data = cellfun(@(x) str2double(strsplit(x, ',')), lines(ind2), 'UniformOutput', false);
    ebar_data = cell2mat(ebar_data);
    ebar_data = ebar_data(:, 2:end-1);
    ebar_data = ebar_data';
    
    time_data = cellfun(@(x) str2double(strsplit(x, ',')), lines(ind+2:end-1), 'UniformOutput', false);
    time_data = cell2mat(time_data);
    time_data = time_data(:, 2:end-1); % cut out integer and NaN
    time_data(:, 1:10) = time_data(:, 1:10) ./ 1.0e3;

    % Covnert time to seconds
    data.t_all = time_data(:, 1);
    data.t_lum = time_data(:, 2);
    data.t_jac = time_data(:, 3);
    data.t_sol = time_data(:, 4);
    data.t_rhs = time_data(:, 5);
    data.t_res = time_data(:, 6);
    data.N_lum = time_data(:, 11);
    data.N_jac = time_data(:, 12);
    data.N_sol = time_data(:, 13);
    data.mean_energy = ebar_data;

end


function c = pair_data(a, b)

    fnames = fieldnames(a);
    for i = 1:numel(fnames)
        assert(isfield(b, fnames{i}), 'b must contain fields from a');
        avar = a.(fnames{i});
        bvar = b.(fnames{i});
        if strcmp('mean_energy', fnames{i})
            c.(fnames{i}) = cat(3, avar, bvar);
        else
            c.(fnames{i}) = [avar, bvar];
        end
    end

end

