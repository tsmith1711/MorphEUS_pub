%% script to implement MRMR
% INPUT: zero_mean_table
%
% OUTPUT: starting_vars

% uses helper mrmr functions 
addpath('./mRMR_0.9_compiled')
addpath('./mRMR_0.9_compiled/mi')


% leave flexible for if using pca_whitened_table or zero_mean_table
chosen_table = zero_mean_table;

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,chosen_table,'OutputFormat', 'uniform');

%find column names
numeric_final_col_names = chosen_table.Properties.VariableNames(numeric_final_data_cols);

% get full numeric values table
numeric_table = chosen_table(:,numeric_final_data_cols);

% get non-numeric values
non_numeric = chosen_table(:,~numeric_final_data_cols);

% get first column (drug column)
drug_list = chosen_table{:,1};

% get the numbers of each unique drug
[~,~,classes] = unique(drug_list);

% get the numerical data
ndata = chosen_table{:,numeric_final_data_cols};

%% discretize data
if disc
    ddata = ndata;

    num_obs = size(ndata,1);

    num_bins = floor(num_obs/5);

    for i = 1:size(ndata,2)
        ddata(:,i) = discretize(ndata(:,i),num_bins);
    end 
else 
    ddata = ndata; 
end 
%% run analysis

% run analysis
kept_cols = mrmr_mid_d(ddata, classes, numvar);

% which columns are kept?
kept_vars = numeric_final_col_names(kept_cols);

% apply this to the data table
reduced_variable_table = horzcat(non_numeric,numeric_table(:,kept_cols)); 

% keep consistent with rest of naming scheme
starting_vars = kept_vars;

% set boolean after_feature_selection to true
after_feature_selection = true; 