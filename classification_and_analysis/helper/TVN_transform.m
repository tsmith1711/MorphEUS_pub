%%% TVN TRANSFORM
% option to include DMSO

% boolean 
if ~exist('after_feature_selection','var')
    after_feature_selection = true;
end 

% if include_DMSO wasn't declared in parent script, set it to false
if ~exist('include_DMSO','var')
   include_DMSO = false;  
end

%% Step 0: Prep the data and get untreated table
% want to normalize the table including both the negative controls and the
% treatments, so first step will get the final_data_table and normalize it

normalized_table = final_data_table; 

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,normalized_table,'OutputFormat', 'uniform');

% get list of drugs and their indicies
drugs = unique(final_data_table.DRUG)';

drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

for i = drugs
    drug = i{1};
    drug_indexes.(drug) = find(strcmp(final_data_table.DRUG, drug)).';
end

% if we have a table where feature selection has already occured 
if exist('starting_vars','var')&& after_feature_selection 
    % get non numeric data
    non_numeric = normalized_table(:,~numeric_final_data_cols); 

    % get the variables that are kept
    kept_variables = normalized_table(:,starting_vars);

    normalized_table = horzcat(non_numeric,kept_variables);

    %find numeric columns of final table
    numeric_final_data_cols = varfun(@isnumeric,normalized_table,'OutputFormat', 'uniform');
end

%find column names
numeric_final_col_names = normalized_table.Properties.VariableNames(numeric_final_data_cols);

%Final check and replace for NaN 
normalized_table = replace_nan(normalized_table, numeric_final_data_cols);

% normalize data
normalized_table{:, numeric_final_data_cols} = normalize(normalized_table{:, numeric_final_data_cols},1,'range'); % normalize by dividing by largest value 

% get values for just untreated controls 
untreated_table = normalized_table(drug_indexes.Untreated,:);

if include_DMSO
    % get values for just DMSO controls
    DMSO_table = normalized_table(drug_indexes.DMSO,:);
end 

%% Step 0.5: prep the data for PCA with just untreated table

if include_DMSO
    % combine untreated and DMSO
    control_table = vertcat(untreated_table,DMSO_table);
else 
    % just use untreated
    control_table = untreated_table;
end 
% get means for each column to use later on full data
col_means = mean(control_table{:, numeric_final_data_cols});

% 0 mean data
control_table{:, numeric_final_data_cols} = bsxfun(@minus,control_table{:, numeric_final_data_cols},mean(control_table{:, numeric_final_data_cols}));

controls = control_table{:, numeric_final_data_cols};


%% messing around with svd

[U,S,V] = svd(cov(controls));

% U is the same as coeff1
% therefore untreated*U = pca_axes
% U*S*V' = cov(untreated)
% want to take out the first couple eigenvalues of S
% to get the ith eigenvalues: S(i,i)

num_pc_removed = 6;

for i = 1:length(num_pc_removed)
    S(i,i) = 0;
    
end 

new_cov = U*S*V';

[U2,S2,V2] = svd(new_cov);

new_pca_axes = controls*U2;

%% Step 1: do PCA on untreated controls

starting_data = control_table{:,numeric_final_data_cols};

[coeff1,pca_axes,vari1] = pca(control_table{:,numeric_final_data_cols});

%% Step 2: Normalize/whiten data
% per-dimension normalization to zero-center and unit variance.
% this is a step I'd like to double check, seemed like he had questions on
% our interpretation of it


pca_axes = (pca_axes - repmat(mean(pca_axes), size(pca_axes,1), 1)) ./ repmat(std(pca_axes), size(pca_axes,1), 1);
axes_size = size(pca_axes);

new_pca_axes = (new_pca_axes - repmat(mean(new_pca_axes), size(new_pca_axes,1), 1)) ./ repmat(std(new_pca_axes), size(new_pca_axes,1), 1);


%%

no_pc1 = pca_axes(:,2:end);

% get the transform matrix from the starting point to the whitented point 
pca_and_whiten_transform =  mldivide(starting_data,pca_axes);
no_pc1_transform = mldivide(starting_data,no_pc1);
new_pca_axes_transform = mldivide(starting_data,new_pca_axes);

%% Step 2.5: prepare full set of data to be transformed
% sort of a continuation of Step 0 

% pseudocode from Mike Ando
% centered_data = subtract_negative_control_mean(original_data)
% whiten_data = pca_and_whiten(centered_data)

% We want to take the means that were subtracted from the untreated
% table to give the columns a mean of zero and subtract them from the full
% set of data

% transform the row of means to be the size of the normalized_table with
% the same row repeating
table_of_means = repelem(col_means,size(normalized_table,1),1);

% subtract the means created with the untreated table from the whole table
normalized_table{:,numeric_final_data_cols} = normalized_table{:,numeric_final_data_cols} - table_of_means;

pca_whitened_table = normalized_table;

% apply pca and whiten transform -- this is what is currently used as the
% final data for analysis
pca_whitened_data = normalized_table{:,numeric_final_data_cols}*pca_and_whiten_transform; 

no_pc1_data = normalized_table{:,numeric_final_data_cols}*no_pc1_transform;

new_pca_axes_data = normalized_table{:,numeric_final_data_cols}*new_pca_axes_transform;

pca_whitened_table{:,numeric_final_data_cols} = pca_whitened_data;

% set after transform to true
%after_TVN = true; 