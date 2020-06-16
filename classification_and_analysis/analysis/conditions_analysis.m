%% Flattens condition data and performs analysis 

%% Init any extra vars
experiments = unique(final_data_table.EXP).';

%% Run the flattening script
flatten_conditions

%% Run analysis
final_data_table = final_combined_table;

%% Sort final_data_table by MoA
% sort by MoA so that drugs with similiar MoA are listed near each other
final_data_table = sort_by_column_name(final_data_table, 'ID');

%% Normalize and initialize table metadata
merged_table_setup

%% Initialize choice_variable and choice_extension
conditions = false;
choice_variable = drugs;
choice_extension = 'DRUG';

%% Run PCA Analysis
PCA_analysis

%% Run KNN analysis
KNN_analysis