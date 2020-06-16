%% Dose response analysis without prompts
%
% DESCRIPTION: Calculates some metadata from prerequisites, runs
% flatten_dose_response to flatten table, and runs PCA and KNN analysis.
%
% PREREQUISITES:
%  - Output of load_workspaces
%  - chosen_ids
%  - Output of get_dose_list
%       - doses_by_id
%       - drugs_by_id
%  - chosen_doses
%  - All prerequisites for PCA and KNN analysis scripts
%
% OUTPUT:
%  - final_data_table flattened across doses
%  - DRUG col is now identical to ID col, and there are a few more cols
%        that describe each row
%
% DEPENDENCIES:
%

%% Filter drugs_by_id and doses_by_id by chosen_doses

%make copy from which to remove
chosen_drugs_by_id = drugs_by_id;
chosen_doses_by_id = doses_by_id;

% Remove drugs/doses that were not chosen
for i=chosen_ids
    id = i{1};
    
    %remove unchosen doses
    if ~isempty(chosen_doses_by_id.(id)) %If empty, don't bother
        not_chosen = ~ismember(chosen_doses_by_id.(id), chosen_doses);
        chosen_doses_by_id.(id)(not_chosen) = [];
    end
    
    %remove unchosen drug-doses
    if length(chosen_drugs_by_id.(id)) > 1 %Don't remove single-dose drugs
        not_chosen = ~contains(chosen_drugs_by_id.(id), chosen_doses);
        chosen_drugs_by_id.(id)(not_chosen) = [];
    end
end
%% Add drugs_to_apply here if desired
%apply = false;

%% Get dose by drug (empty entry for single dose drugs)
chosen_doses_by_drug = struct();
for i=chosen_ids
    id = i{1};
    
    chosen_doses_by_drug.(id) = struct();
    
    for j=1:length(chosen_drugs_by_id.(id))
        drug = chosen_drugs_by_id.(id){j};
        
        if isempty(chosen_doses_by_id.(id))
            chosen_doses_by_drug.(id).(drug) = [];
        else
            chosen_doses_by_drug.(id).(drug) = chosen_doses_by_id.(id){j};
        end
    end
end
%% Flatten final_data_table

flatten_dose_response

%% Run analysis
final_data_table = final_combined_table;

%% Sort final_data_table by MoA
% sort by MoA so that drugs with similiar MoA are listed near each other
final_data_table = sort_by_column_name(final_data_table, 'ID');


%% Regenerate DRUG_EXP

%%% reorganize table by numeric/non numeric %%%
%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');
%find column names
numeric_final_col_names = final_data_table.Properties.VariableNames(numeric_final_data_cols);
% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);
% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols); 
%%%

% create drug-experiment column
drug_exp = table();
drug_exp.DRUG_EXP = final_data_table.ID_EXP;

% reorganize table with new column
final_data_table = [non_numeric drug_exp numeric_data];

%% initialize table metadata
merged_table_setup

%% Initialize choice_variable and choice_extension
% If multiple workspaces, drugs will be denoted by drug_exp
% If only one workspace, drugs will be denoted by drug

if conditions %if multiple workspaces
    % get list of drug-experiments from table, to use instead of "drugs"
    drug_exps = unique(final_data_table.DRUG_EXP,'stable')';
end 

% if choosing more than one workspace, need to use drug-experiment instead
% of drugs for choice variable
if conditions
    choice_variable = drug_exps;
    choice_extension = 'DRUG_EXP';
else
    choice_variable = drugs;
    choice_extension = 'DRUG';
end
%     
% 
% %% Run PCA Analysis
%chosen_drugs = chosen_ids;

% PCA_analysis
% 
% %% Run KNN analysis
% 
% KNN_analysis