function output_table = normalize_and_zero_mean(table,cols)
%NORMALIZE_AND_ZERO_MEAN Executes a multi-step perparation on the
%provided table across each of the specified columns.
%   INPUT: table - the table which is to be normalized
%          cols  - a horizontal logical vector of same length as the number
%          of features in the table which specifies which features should
%          be normalized
%   DESCRIPTION:
%          For each column specified, Replaces NaN's with 0, divides by the
%          maximum value in the column to normalize between 0 and 1, and
%          zero means.
%   OUTPUT: output_table - the normalized table.

%% Check for bad input

%Check if same length
if size(cols, 2) ~= width(table)
    error('Input Error: table and cols must be of same size.')
end

%Check that is logical (ones and zeros only)
if sum(cols ~= 1 & cols ~= 0) > 0
    error('Input Error: cols must be a logical vector containing only 1 or 0')
end

%make sure all of cols is numeric
numeric_cols = varfun(@isnumeric, table,'OutputFormat', 'uniform');
if sum(~numeric_cols & cols) > 0
    error('Input Error: cols must specify only numeric columns in table')
end
    
%% Perform Normalization

output_table = table;

%Replace the NaN's with 0
output_table = replace_nan(output_table, cols);
% normalize data
output_table{:, cols} = normalize(output_table{:, cols},1,'range'); % normalize by dividing by largest value 
% 0 mean data
output_table{:, cols} = bsxfun(@minus, output_table{:, cols}, mean(output_table{:, cols}));

end

