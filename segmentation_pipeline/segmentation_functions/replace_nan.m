function out_table = replace_nan(table, numeric_cols)
    if(size(table,1) == 0)
        out_table = table;
    else
        labels = table(:, ~numeric_cols);
        data = table(:, numeric_cols);
        no_nan_data = fillmissing(data, 'constant', 0);
        out_table = [labels no_nan_data];
    end
end