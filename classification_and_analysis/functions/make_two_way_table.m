function two_way_table = make_two_way_table(data_table,colvar,rowvar)
% Input a table and two variables that you would like to make into a table
% 

    all_colvars = unique(data_table.(colvar),'stable')';
    
    all_rowvars = unique(data_table.(rowvar))';

    if length(all_colvars) == length(all_rowvars) || length(all_rowvars) == length(all_colvars)-1 || length(all_rowvars)-1 == length(all_colvars)
         all_rowvars = all_colvars;
        
    end 
    empty_array = zeros(length(all_rowvars),length(all_colvars));
    
    two_way_table = array2table(empty_array,'VariableNames',all_colvars,...
        'RowNames',all_rowvars);
    
    
    for i = all_colvars 
        col = i{1};

        row_inds = find(strcmp(data_table.(colvar),col));

        drug_rows = data_table(row_inds,:);

        neighbor_list = drug_rows.(rowvar);
        for j = all_rowvars
            row = j{1};
            if ismember(row,neighbor_list)
                % find where it is in the drug rows
                n_ind = find(strcmp(drug_rows.(rowvar),row));
                weight = drug_rows.DISTANCE(n_ind);
                two_way_table{row,col} = weight;
            else
                two_way_table{row,col} = 0;
            end 

        end 


    end 

end

