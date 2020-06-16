function [drug_medians,drugs] = table_median(knn_table,choice_extension)
% Creates a new table that takes the median values for the last table
% helper for knn_helper 
% knn_table: whatever your table of choice that you want to get the median
% values for and perform knn
% choice_extension: for when you select either one worksapce or more than
% one workspace, will be DRUG or DRUG_EXP depending 
% outputs new drug_medians table and the list of drugs in the table

    knn_drugs = knn_table.(choice_extension);
    
    numeric_final_data_cols = varfun(@isnumeric,knn_table,'OutputFormat', 'uniform');
    
    knn_data = knn_table{:,numeric_final_data_cols};
    drugs = unique(knn_drugs,'stable')';
    
    choice_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);
    for i = drugs
        drug = i{1};
        choice_indexes.(drug) = find(strcmp(knn_table.(choice_extension), drug)).';
    end
    
    % preallocate data: # rows of drugs, # cols of knn_data
    drug_medians = zeros(length(drugs),size(knn_data,2));
    % initialize counter 
    ct = 1;
    for i = drugs
        % get individual drug
        drug = i{1};

        % new row for each median -- will be in same order as drugs veector
        drug_medians(ct,:) = median(knn_data(choice_indexes.(drug),:));

        % increasee row countere 
        ct = ct+1;
    end

end

