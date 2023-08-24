%% Function for reading cross sections
function proc = read_lxcat_file(filename)
    
    % Open file and read to cell array
    file_lines = extract_lines_from_file(filename);
    
    % Cut the file into cell array of process blocks
    block_text = cut_into_raw_text_blocks(file_lines);
    
    % Loop and fill
    i = 1;
    proc(i) = parse_process(block_text{i,1}, block_text{i,2}); 
    proc = repmat(proc(i), size(block_text, 1), 1);
    for i = 2:size(block_text, 1)
        proc(i) = parse_process(block_text{i,1}, block_text{i,2}); 
    end

    % Add duplicates for superelastic processes
    proc = add_superelastics(proc);
    
 end


 %% Create new processes for super-elastics
 function proc = add_superelastics(proc)
    
    % Locate reversible Processes
    is_rev = find([proc(:).is_reversible]);

    % Copy those and append to end of proc stack
    Ns = numel(proc);
    proc = [proc; proc(is_rev)];
    
    % Go through and correct details in the 
    for i = (Ns+1):numel(proc)

        % Inidicate this process is the superelastic inverse
        proc(i).is_superelastic = true;
        
        % Locate Linked process 
        j = is_rev(i-Ns);
        proc(i).linked_process = j;
        proc(j).linked_process = i;
        
        % Swap products and reactants
        proc(i).species_name  = proc(j).product_names;
        proc(i).product_names = proc(j).species_name;
        
        % Generate Super-elastic cross section data
        %   New data is tabulated at eps-de
        eps_j = proc(i).eps_j - proc(i).de;
        sigma_j = proc(i).sigma_j .* proc(i).eps_j ./ (proc(i).gratio .* eps_j);
        
        % Clear any negatives (from super-elastic shift)
        is_oob = eps_j <= 0.0;
        eps_j(is_oob) = []; 
        sigma_j(is_oob) = [];
        
        % Overwrite interpolation functions
        [proc(i).sigma_linear, proc(i).sigma_log] = store_interp_func(eps_j, sigma_j);
        
        % Store the new eps_j, sigma_j
        proc(i).eps_j = eps_j;
        proc(i).sigma_j = sigma_j;

    end

 end


%% Open file and read lines
function file_lines = extract_lines_from_file(filename)

    % Open File
    fid = fopen(filename);

    % Read file into cell array of strings
    file_lines = textscan(fid, "%s", "Delimiter", "\n");
    file_lines = file_lines{1};
    fclose(fid);
    
end


%% Parses file into blocks for each process,
%   using --- lines and header strings as edges. 
%   returns block_text cell array of Ns x 2 size, first column is header, 
%   second column is values between *** lines
function block_text = cut_into_raw_text_blocks(file_lines)
    
    % Locate dash lines
    ind_dash = reshape( find(strncmp(file_lines, '---', 3)), [], 1); % as row vector
    
    % Check there is an even number of ---
    assert(~mod(numel(ind_dash),2), 'odd number of --- lines detected in file');
    
    % Reshape into Ns x 2
    ind_dash = reshape(ind_dash, 2, [])';

    % Species count, init ptext cell array
    Ns = size(ind_dash, 1);
    block_text = cell(Ns, 2);

    % Locate Process Strings to mark start
    % This assumes headers are all caps, first non-whitespace entry in a
    % given line
    process_type_str = {'MOMENTUM', 'ELASTIC', 'EXCITATION', 'IONIZATION', 'ATTACHMENT'};
    block_start = find(startsWith(strtrim(file_lines), process_type_str)); 
    assert(numel(block_start)==Ns, 'Nepsmber of header strings must match number of dash lines divided by 2')
    
    % Loop process and seperate
    for i = 1:Ns
        block_text(i, 1) = {file_lines(block_start(i):ind_dash(i, 1)-1)};
        block_text(i, 2) = {file_lines(ind_dash(i, 1)+1:ind_dash(i, 2)-1)};
    end

end


%% Process block of text
function proc = parse_process(header_txt, data_txt)
    
    %% Parse Type
    process_type_str = {'EFFECTIVE', 'ELASTIC', 'EXCITATION', ...
                        'IONIZATION', 'ATTACHMENT'};
    proc.collision_type_str = strtrim(header_txt{1});
    cid = find(strcmp(process_type_str, proc.collision_type_str ), 1);
    
    % Useful logic (if redundant)
    proc.collision_type = cid;
    proc.is_effective = cid == 1;
    proc.is_elastic = cid == 2;
    proc.is_excitation = cid == 3;
    proc.is_ionization = cid == 4;
    proc.is_attachment = cid == 5;
    
    % Other Useful logic
    proc.is_elastic_or_effective = cid <=2;
    proc.is_inelastic = cid > 2;
    proc.is_nonconservative = cid > 3;
    proc.is_source_type = cid ==3 || cid==4;
    
    %% Default values that may not be filled
    proc.mratio = NaN;
    proc.de = NaN;
    proc.is_reversible = false;
    proc.is_superelastic = false;
    proc.linked_process = 0;
    proc.gratio = NaN;
    
    %% Parse Header lines 2-3 (reaction name and de/gratio/mratio)
    if proc.is_elastic_or_effective

        % Effective or Elastic
        proc.reaction_name = parse_comment('PROCESS', header_txt);
        proc.species_name = strtrim(header_txt{2});
        proc.mratio = str2double(header_txt{3});
        proc.product_names = {proc.species_name};
        
    else
        
        % Get Reaction Name string
        proc.reaction_name = header_txt{2};
        
        % Ensure it has a reaction symbol -> or <->
        assert(contains(proc.reaction_name, '->'), 'Bad read of reaction string')
        proc.is_reversible = contains(proc.reaction_name, '<->');
        
        % Parse Reactants
        full_names = strtrim(split(proc.reaction_name, {'->','<->'}));
        reac_names = strtrim(split(full_names{1}, ' + ')); 
        assert(numel(reac_names)==1, 'More than one reactant species detected!');
        proc.species_name = reac_names;
        
        % Parse Product Species
        proc.product_names = strtrim(split(full_names{2}, ' + ')); 
        if proc.is_reversible
            assert(numel(proc.product_names)==1, 'Reversible process cannot have multiple products!')
        end
        
        % If not attachment, get threshold energy and ratio of degeneracies
        if ~proc.is_attachment
            
            % Split and covnert line to double
            values = str2double(split(strtrim(header_txt{3})));
            
            % First value is assumed to be threshold in eV
            proc.de = values(1);
            
            % Second value is g2/g1, default is 1
            if numel(values) > 1
                proc.gratio = values(2);
            else % Default g2/g1 = 1
                proc.gratio = 1.0;
            end

        end
        
    end
    
    %% Parse Comment Block
    proc.comment_species = parse_comment('SPECIES', header_txt);
    proc.comment_process = parse_comment('PROCESS', header_txt);
    proc.comment_param   = parse_comment('PARAM',   header_txt);
    proc.comment_other   = parse_comment('COMMENT', header_txt);
    proc.comment_columns = parse_comment('COLUMNS', header_txt);
    
    % Parse Dates
    proc.date_updated = parse_comment('UPDATED', header_txt);
    proc.date_used = datetime('now');
    
    %% Process Data
    
    % Extract data from data block
    data = reshape(sscanf(sprintf('%s ',data_txt{:}),'%f',[1, Inf]), 2, [])';
    proc.eps_j = data(:,1);
    proc.sigma_j = data(:, 2);
    
    % Store interpolation functions
    [proc.sigma_linear, proc.sigma_log] = store_interp_func(proc.eps_j, proc.sigma_j);
    
end


%% Interpolation Functions
function [sigma_linear, sigma_log] = store_interp_func(eps_j, sigma_j)

    % linear Interp Model
    sigma_linear = griddedInterpolant(eps_j, sigma_j, 'linear', 'none');
    
    % Apply log-log
    eps_j(~logical(eps_j)) = 1e-30; % Avoid nans!
    eps_j = log10(eps_j);
    sigma_j = log10(sigma_j);
    sigma_log = griddedInterpolant(eps_j, sigma_j, 'linear', 'none');
    
end


%% Parse Comments
function results = parse_comment(keyword, header)
    
    % Clean up keyword
    keyword = strtrim(keyword);
    if ~endsWith(keyword, ':')
        keyword = [keyword, ':'];
    end
    
    % Find in blocks
    ind = startsWith(strtrim(header), keyword);
    
    % Process to single string
    results = strjoin(header(ind));
    results = erase(results, keyword);
    results = strtrim(results);
    
    % Get actual data
    if contains(keyword, 'UPDATED')
        results = datetime(results);
    end

end
