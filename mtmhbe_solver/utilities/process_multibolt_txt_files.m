
function data = process_multibolt_txt_files(varargin)
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
    ind = find(strcmp(lines, ' SUMMARY:' ));
    
    % Find Nsim
    Nsim = str2double( lines{ind+1} );

    ebar_data = cellfun(@(x) str2double(strsplit(x, ',')), lines(1:ind-1), 'UniformOutput', false);
    ebar_data = cell2mat(ebar_data);
    ebar_data = ebar_data(:, end);
    ebar_data = reshape(ebar_data, 55, []);

    Nloop = size(ebar_data, 2);

    ind_offset = 2;
    
    time_data = cellfun(@(x) str2double(strsplit(x, ',')), lines(ind+1+ind_offset:end), 'UniformOutput', false);
    time_data = cell2mat(time_data);
    time_data = time_data(:, 2:end);
    
    
    % Patch Timers
    time_data(:, 3:6) = time_data(:, 3:6) ./ (Nsim * 1e6);

    % Peal Out Times
    data.Nsolve = time_data(:, 1);
    data.Ngrid = time_data(:, 2);
    data.total = time_data(:, 3);
    data.jac_solve = time_data(:, 6);
    data.jac_solve_per = data.jac_solve ./ data.Nsolve;
    data.jac_grid = time_data(:, 5);
    data.jac_grid_per = data.jac_grid ./ data.Ngrid;
    data.jac   = data.jac_solve + data.jac_grid;
    data.solve = time_data(:, 4);
    data.resid = data.total - data.jac - data.solve;
    data.mean_energy = ebar_data;
    data.Nloop = Nloop;

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