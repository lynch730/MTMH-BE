function b2 = combine_bolsig(b)

    b2 = b(1);
    for i = 2:numel(b)
        fn = fieldnames(b2.settings);
        for j = 1:numel(fn)
            b2.settings.(fn{j})(:, end+1) = b(i).settings.(fn{j})(:, 1);
        end
        fn = fieldnames(b2.moments);
        for j = 1:numel(fn)
            value = b(i).moments.(fn{j})(:, 1);
            if isempty(value)
                value = 0;
            end
            b2.moments.(fn{j})(:, end+1) = value(:, 1);
        end
        fn = fieldnames(b2.rates);
        for j = 1:numel(fn)
            b2.rates.(fn{j})(:, :, end+1) = b(i).rates.(fn{j})(:, :, 1);
        end
        b2.eedf(:, i) = b(:, i).eedf;
    end

%     for i = 1:numel(b(1).species_names)
%         if numel(b(1).species_fractions{i})>1
%             b2.variables.(['spec_frac_', b.species_names{i}]) = b.species_fractions{i};
%         end
%     end

end