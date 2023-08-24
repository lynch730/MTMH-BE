function mdata = extract_multibolt(mpath, EN_array)
    
    N_EN = numel(EN_array);

    % Get list of all files in root dir
    
    %% Extract Gas Data
    spec_details =  'species_indices.txt';
    mdata.PerGas.species = read_file_species(fullfile(mpath, spec_details));
    
    for i = 1:numel(mdata.PerGas.species.species_idx)
        idx = mdata.PerGas.species.species_idx(i);
        spec_folder = fullfile(mpath, 'PerGas', ['Gas_',num2str(idx, '%i')]);
        mdata.PerGas.species.species_idx(i) = idx + 1;
        sfiles = dir(spec_folder);
        sname = ['spec_', num2str(i, '%i')];
        mdata.PerGas.(sname).k_ela = zeros(N_EN, 1);
        mdata.PerGas.(sname).k_exc = zeros(N_EN, 1);
        mdata.PerGas.(sname).k_iz = zeros(N_EN, 1);
        mdata.PerGas.(sname).k_att = zeros(N_EN, 1);
        for j = 1:numel(sfiles)
            if sfiles(j).isdir
                continue
            end
            str = split(sfiles(j).name, {'_', '.'});
            kcnt = str2double(str{end-1}) + 1;
            ktable = read_file(fullfile(spec_folder, sfiles(j).name));
            switch str{end-3}
                case 'ela'
                    mdata.PerGas.(sname).k_ela(:, kcnt) = ktable.k_ela_N;
                case 'exc'
                    mdata.PerGas.(sname).k_exc(:, kcnt) = ktable.k_exc_N;
                case 'iz'
                    mdata.PerGas.(sname).k_iz(:, kcnt) = ktable.k_iz_N;
                case 'att'
                    mdata.PerGas.(sname).k_att(:, kcnt) = ktable.k_att_N;
                case 'sup'
                    mdata.PerGas.(sname).k_sup(:, kcnt) = ktable.k_sup_N;
                otherwise 
                    continue
            end
            
        end
    end
    
    %% Total rates
    Tpath = fullfile(mpath, 'Total');
    Tfiles = dir(Tpath);
    for j = 1:numel(Tfiles)
        if Tfiles(j).isdir
            continue
        end
        str = split(Tfiles(j).name, {'_', '.'});
        ktable = read_file(fullfile(Tpath, Tfiles(j).name));
        switch str{end-2}
            case 'ela'
                mdata.Total.k_ela = ktable.k_ela_N;
            case 'exc'
                mdata.Total.k_exc = ktable.k_exc_N;
            case 'iz'
                mdata.Total.k_iz = ktable.k_iz_N;
            case 'att'
                mdata.Total.k_att = ktable.k_att_N;
            case 'sup'
                mdata.Total.k_sup = ktable.k_sup_N;
            otherwise 
                continue
        end
    end

    %%  EEDF F0
    Epath = fullfile(mpath, 'EEDFs_f0');
    Efiles = dir(Epath);
    for j = 1:numel(Efiles)
        if Efiles(j).isdir 
           continue
        end
        ktable = read_file(fullfile(Epath, Efiles(j).name));
        str = split(Efiles(j).name, {'_', '.'});
        kcnt = str2double(str{end-1}) + 1;
        mdata.F0.eV(:, kcnt) = ktable.eV;
        mdata.F0.f0(:, kcnt) = ktable.f0;
    end

    %%  EEDF F1
    Epath = fullfile(mpath, 'EEDFs_f1');
    Efiles = dir(Epath);
    for j = 1:numel(Efiles)
        if Efiles(j).isdir 
           continue
        end
        ktable = read_file(fullfile(Epath, Efiles(j).name));
        str = split(Efiles(j).name, {'_', '.'});
        kcnt = str2double(str{end-1}) + 1;
        mdata.F1.eV(:, kcnt) = ktable.eV;
        mdata.F1.f0(:, kcnt) = ktable.f1;
    end
    
    %% 
    files = dir(mpath);
    excluded_files = {'sweep_indices.txt', spec_details, 'RUN_DETAILS.txt'};
    for i = 1:numel(files)
        if files(i).isdir || contains(files(i).name, excluded_files)
            continue
        end
        
        ktable = read_file(fullfile(mpath, files(i).name));
        fname = erase(files(i).name, '.txt');
        mdata.(fname) = ktable{:, 2};

    end
    
    mdata.EN_array = EN_array;
    mdata.N_EN = numel(EN_array);
    
end

function table = read_file_species(filename)
    table = read_file(filename);
    if numel(table.Properties.VariableNames) > 2
        table2.species_idx = table.Var1;
        for i = 1:numel(table.Var1)
            str = [];
            for j = 2:numel(table.Properties.VariableNames)
               value = table.(['Var', num2str(j, '%i')]);
               str = [str, ' ', value{i}];
            end
        	table2.species_name{i} = str;
        end
        table = table2;
    end
end

function table = read_file(filename)
    table = readtable(filename , 'CommentStyle',{'#'});
end

