%% Helper Script to renormalize final_data_table and redeclare table metadata in a merged workspace
% This script needs to be called when concatenating multiple workspaces together in order to
% renormalize to the new concatenated space and recalculate the metadata
% variables which describe the table.
%
% Called when using multiple workspaces in load_workspaces.m
%
% DESCRIPTION: - *NO NOT ANYMORE* Calls normalize_and_zero_mean() on final_data table. 
%              - Initializes several variables which describe the data in
%                the new normalized table. Many of these exist beforehand,
%                we just need to make sure they reflect the concatenated
%                and "prepared" data
%
% INPUT:  - final_data_table
%
% OUTPUT: - final_data_table, now normalized and zero meaned
%         - "Metadata variables" which describe the now normalized table
%                 - numeric_final_data_cols
%                 - numeric_fianl_col_names
%                 - normalized_table
%                 - ndata: a matrix of numeric, 
%                 - drugs
%                 - drug_indexes
%                 - imgs
%                 - img_indexes
%                 - yvalues: drug labeles for clustergram

%% Find column information

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');
%find column names
numeric_final_col_names = final_data_table.Properties.VariableNames(numeric_final_data_cols);

%% This section has been removed
% We don't want to overwrite final_data_table
% Normalization is still always done for us in PCA_analysis

% %% Normalize table
% 
% final_data_table = normalize_and_zero_mean(final_data_table, numeric_final_data_cols);
% 
% %% Get the numeric data from the table
% 
% % get numeric data from normalized table as a matrix for later use in PCA calculations
% ndata = final_data_table{:,numeric_final_data_cols}; 

%% getting unique non-numeric values and indicies, and initializing a few more variables

drugs = unique(final_data_table.DRUG)';

%find row indexes for each drug
drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);
for i = drugs
    drug = i{1};
    drug_indexes.(drug) = find(strcmp(final_data_table.DRUG, drug)).';
end

% get unique images
imgs = unique(final_data_table.IMG)';

% get row indexes for each image %% COMMENTED OUT FOR DOSE RESPONSE %%
% img_indexes = cell2struct(cell(1,length(imgs)), imgs, 2);
% for i = imgs
%     drug = i{1};
%     img_indexes.(drug) = find(strcmp(final_data_table.IMG, drug)).';
% end



% get a numerical list for each drug (eg if there are 3 Cer rows, will
% create a cell array including values Cer1 Cer2 Cer3). This is used in the
% creation of a clustergram. 
yvalues = (final_data_table.DRUG);
for i = drugs
    drug = i{1};
    for j = 1:length(drug_indexes.(drug))
        indx = drug_indexes.(drug)(j);
        yvalues(indx) = strcat(yvalues(indx),num2str(j));
    end
end

% get non-numeric information
non_numeric = final_data_table(:,~numeric_final_data_cols);
