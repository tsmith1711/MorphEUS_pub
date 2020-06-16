function fixed_big_data = column_validation(all_spaces,prompt_user)
% during the step that was big_data = struct2cell(all_spaces)
% validate that the tables in the given structure all have the same columns
% and, if not, prompt user to delete some in order to merge
% fixed_big_data is a cell array created from the structure of tables
% prompt is variable to prompt user or if its okay to remove variables,
% false means variables will be auto removed
% default is false

    if ~exist('prompt_user','var')
        prompt_user = false;
    end 


    tables = fieldnames(all_spaces)';
    % create a struct to store all column names from all tables
    table_cols = struct;
    % go through each table in structure
    ct = 1;
    for i = tables
        current_field = i{1};
        current_table = all_spaces.(current_field);
        % store the variable names 
        table_cols.(current_field) = current_table.Properties.VariableNames;
        ct = ct+1;
    end 
    
    % check if all values in # column names are the same. If they are, we
    % are good to go, if they are not, we gotta fix some things. 
    
    % if all the values are not equal to the first value 

    % time for some fixing 

    % setdiff needs to use the largest one first 
    % we need to loop through two at a time... convert to cell array? 
    cell_array =  struct2cell(table_cols);
    diff_vars = [];
    for i = 1:length(cell_array)-1
        A = cell_array{i};
        B = cell_array{i+1};
        % order for setdiff matters so do both 
        diff_vars_1 = setdiff(A,B);
        diff_vars_2 = setdiff(B,A);

        % on first time through, just add to the diff_vars
        if isempty(diff_vars)
            diff_vars = unique([diff_vars_1 diff_vars_2]);
        % if following times through, want to make sure to keep the
        % diff_vars detected last time 
        else 
            new_diff_vars = unique([diff_vars_1 diff_vars_2]);
            diff_vars = unique([diff_vars new_diff_vars]);
        end

    end 

    %%
    % now that we've gotten all of the differeneces, need to reemove
    % these differences 
    if~isempty(diff_vars)
        rmv_str = "";
        for j = diff_vars
            drug = j{1};
            rmv_str = strcat(rmv_str, drug, " ");
        end
        if prompt_user
            rmv = questdlg(strcat("remove variables", rmv_str, "?"));
        else
           rmv = 'Yes'; 
        end
        switch rmv
            case 'Yes'
                for i = tables
                    current_field = i{1};
                    % removevars will throw an error if that variable is
                    % already removed, so loop through one at a time with a
                    % try catch I guess
                    for w = diff_vars
                        try
                            all_spaces.(current_field) = removevars(all_spaces.(current_field),w);
                            disp(strcat("removed ", w{1}, " from ", current_field))
                        catch
                            disp(strcat(w{1}, " was not in ", current_field))
                        end 
                    end 

                end 

            case 'No'
                error(strcat("Could not concatenate tables due to ", rmv_str))
        end
    end 
    fixed_big_data = struct2cell(all_spaces); 
end

