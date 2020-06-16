function sorted_table = sort_by_column_name(table, colname)
%SORT_BY_COLUMN_NAME Sorts rows of table by column name
%   Calls sortrows on the column name

colnames = table.Properties.VariableNames;
col_index = find(strcmp(colname, colnames));

sorted_table = sortrows(table, col_index);
end

