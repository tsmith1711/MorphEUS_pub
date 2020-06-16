function merged_table = mergeImages(struct, fields)
    filtered_fields = fields(isfield(struct, fields));
    tables_to_merge = cellfun(@(f)struct.(f),filtered_fields,'UniformOutput',false);
    
    %remove rownames to prevent conflict
    for i = 1:size(tables_to_merge, 2)
        tables_to_merge{i}.Properties.RowNames = {};
    end
    
    merged_table = vertcat(tables_to_merge{:});
end