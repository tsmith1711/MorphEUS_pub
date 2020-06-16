function final_data_table = combine_workspaces(chosen_workspaces,workspace_directory)
% helper function for load_workspaces when doing more than one workspace 
% input chosen workspaces and workspace_directory, get the combined
% final_data_table 
    % create a structure to store the final_data_table s from each workspace
    all_spaces = struct;
    for i = chosen_workspaces
       % get each individual spacee
       one_space = i{1};
       % need to get final data table from each space and concat final_data_table
       fpath = strcat(workspace_directory, '/', one_space);
       % get the name without the .mat extension
       before_mat = extractBefore(one_space,".");
       
       %Make sure doesnt have the variables we want to use
       if any(strcmp(who('-file', fpath), 'all_spaces')) || any(strcmp(who('-file', fpath), 'big_data'))
           error(strcat("Vars that will break load_workspaces exist in ", before_mat))
       end
       
       % load the final_data_table from the workspace
       load(fpath) % Might be able to load just final_data_table, after_TVN
       %save the final_data_table under the name of the workspace
       all_spaces.(before_mat) = final_data_table;
    end

    % convert into one large final_data_table
    big_data = column_validation(all_spaces);

    final_data_table = vertcat(big_data{:});
    
end

