%% Helper Script for PCA_analysis.m
%%% Validates workspace variables chosen_drugs and drugs_to_apply, and
%%% filters final_data_table so that only the chosen drugs are included.
%%%
%%% INPUT:
%%% Intended to only be called internally by PCA_Analysis.m, but in theory
%%% the only prerequisites are:
%%%     - chosen_drugs
%%%     - drugs_to_apply
%%%     - the output of load_workspaces.m (final_data_table + metadata vars)
%%%
%%% OUTPUT: final_data_table, with the chosen drugs removed
%%%         apply , whether there are drugs to apply
%%%         drug_to_apply, a subsection of the final_data_table only for the drug to apply 
%%%         choice_indexes
%%%         chosen chosen_drugs
%%%         

%% Verify and initialize drug selections
% set default to false for this: if the user selects untreated then nothing
% will happen here but it needs a default 

if ~exist('use_medians_only', 'var')
    use_medians_only = false;
    disp('use_medians_only defaulted to false')
end 

% Make sure list for drugs to use exists
if ~exist('chosen_drugs', 'var')
    error('Fatal Error: Must specify chosen_drugs for PCA_analysis.m')
end

% Make sure chosen_drugs are all in the drug list
%%% choice_variable is either drugs or drug_exps, depending on whether
%%% there were multiple workspaces or not

for i = chosen_drugs
    chosen_drug = i{1};
    if ~any(strcmp(choice_variable, chosen_drug))
        error(strcat("Fatal Error: Drug ", chosen_drug, " in chosen_drugs is not contained in drugs list."))
    end
end

% apply = true if drugs_to_apply exists and has elements, apply=false otherwise
if ~exist('drugs_to_apply', 'var') || length(drugs_to_apply) <= 0
    apply = false;
else
    apply = true;
    
    % Also, if we are applying drugs, make sure all of them exist in the drugs list
    for i = drugs_to_apply
    chosen_drug = i{1};
        if ~any(strcmp(choice_variable, chosen_drug))
            error(strcat('Fatal Error: Drug ', chosen_drug, ' in drugs_to_apply is not contained in drugs list.'))
        end
    end
end

% find indexes of each drug in the final data table 
choice_indexes = cell2struct(cell(1,length(choice_variable)), choice_variable, 2);
for i = choice_variable
    drug = i{1};
    choice_indexes.(drug) = find(strcmp(final_data_table.(choice_extension), drug)).';
end

% duplicate for drugs chosen to train PCA
chosen = chosen_drugs;

%%
% want to make them keep untreated, put in a boolean to get rid of it later
% after TVN 
% force user to select untreated
if ~any(contains(chosen,'Untreated'))
    chosen = [chosen {'Untreated'}];
    remove_after_TVN = true;
    disp("in filter_drugs: remove_after_TVN set to true")
end 

if ~exist('remove_after_TVN','var')
    remove_after_TVN = false; 
end 

% if user selected drug(s) to apply 
if apply 
    % Duplicate for drugs to apply
    addedDrug = drugs_to_apply;

    % get indicies of all drugs selected to apply 
    all_kept_inds = [];
    for i = addedDrug
        drug = i{1};
        all_kept_inds = [all_kept_inds choice_indexes.(drug)];
    end 

    % get drug(s) selected to apply
    drug_to_apply = final_data_table(all_kept_inds,:);
end 

% drugs user decided to ignore
removed = choice_variable;
indx = ismember(choice_variable, chosen_drugs);
removed(indx) = [];
clear indx;
% the applied drug should be in removed 

% make sure that Untreated is not removed
if any(contains(removed,'Untreated'))
    untreated_loc = contains(removed,'Untreated');
    removed = removed(~untreated_loc);
    remove_after_TVN = true;
     disp("in filter_drugs: remove_after_TVN set to true again")
end 

% remove ignored drugs from the list 
removallist = [];
for i = removed
    drug = i{1};
    removallist = [choice_indexes.(drug) removallist];
    
end
final_data_table(removallist,:) = [];

% if we want to only use median values and no q1, q3, iqr
if use_medians_only
    % find q1 columns
    q1_cols = contains(final_data_table.Properties.VariableNames,'_q1');
    q1_col_names = final_data_table.Properties.VariableNames(q1_cols);
    % find q3 columns
    q3_cols = contains(final_data_table.Properties.VariableNames,'_q3');
    q3_col_names = final_data_table.Properties.VariableNames(q3_cols);
    % find IQR cols
    iqr_cols = contains(final_data_table.Properties.VariableNames,'_v');
    iqr_col_names = final_data_table.Properties.VariableNames(iqr_cols);
    
    all_cols_to_remove = [q1_col_names q3_col_names iqr_col_names];
    
    final_data_table = removevars(final_data_table,all_cols_to_remove);
end 

numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');