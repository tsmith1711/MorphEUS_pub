function out_table = create_exp_column(data_table, experiment)
%CREATE_EXP_COLUMN Creates an EXP field in a table, and fills with the
%given EXP
%   Inserts EXP column at end of non numeric fields

if ~ischar(experiment)
    error('Experiment name must be character array (single quote string)')
end

% Split current table
numeric_cols = varfun(@isnumeric, data_table,'OutputFormat', 'uniform');
numeric_table = data_table(:, numeric_cols);
non_numeric_table = data_table(:, ~numeric_cols);

% Make exp col
exp_data = cell(size(data_table, 1), 1);
exp_data(:) = {experiment};
exp_col = table();
exp_col.EXP = exp_data;

% Insert col
out_table = [non_numeric_table exp_col numeric_table];
end

