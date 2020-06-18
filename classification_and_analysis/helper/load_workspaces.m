%% Loads and Prepares a list of workspaces for PCA_analysis
% DESCRIPTION:
%   Verifies and loads workspaces in chosen_workspaces and vertically
%   concatenates them.
%
%   Calls merged_table_setup.m if multiple tables were merged (see
%   that script's description)
%
% PREQUISITES:
%   - chosen_workspaces exists in the MATLAB workspace, containing a
%       horizontal cell array of string names for workspaces (including
%       file extension)
%   - (OPTIONAL) remove_INH_control - whether to run remove_INH_control on
%       on final data table (default: true)
%   - (OPTIONAL) remove_extra_controls - Will remove:
%       EtOH, MeOH, NaOH, water (default: true)
%
%   - (OPTIONAL) workspace_directory - the directory where to look for
%   workspaces (default: './workspaces')
%
%
% OUTPUT
%   - Throws error if chosen_workspaces does not exist
%   - Prints the workspaces that are included from chosen_workspaces, and
%       throws an error if any one is found that does not exist in ./workspaces
%  Adds to workspace:
%   - final_data_table, concatenated from all workspaces, and
%     and renormalied/zero-meaned if merged
%   - All the metadata outputs of merged_table_setup.m
%   - Other metadata vars:
%       - conditions - true if multiple workspaces, false if only one
%       - choice_extension - either 'DRUG' or 'DRUG_EXP' depending on var conditions
%       - choice_variable - a list of either drugs or drug_exps depending on var conditions

%% Generate defaults
if ~exist('remove_INH_control', 'var')
    remove_INH_control = true;
    disp("Defaulted remove_INH_control to TRUE in load_workspaces")
end

if ~exist('remove_extra_controls', 'var')
    remove_extra_controls = true;
    disp("Defaulted remove_extra_controls to TRUE in load_workspaces")
end

if ~exist('remove_bad_treated', 'var')
    remove_bad_treated = true;
    disp("Defaulted remove_bad_treated to TRUE in load_workspaces")
end

if ~exist('workspace_directory', 'var')
    workspace_directory = './workspaces';
    disp("Defaulted workspace_directory to './workspaces' in load_workspaces")
end

%% Import paths
addpath(workspace_directory)
%% Get all existing workspaces and check list against them

%Make sure chosen_workspaces exists
if ~exist('chosen_workspaces', 'var') || length(chosen_workspaces)  <= 0
    error('Must define variable chosen_workspaces before running load_workspaces.m')
end

% get all workspaces in the workspace folder
workspace_list = get_workspaces(workspace_directory);

% Print all workspaces to console, and throw error if any one doesn't exist in ./workspaces/
disp('Importing workspaces:')
for i = chosen_workspaces
    wksp = i{1};
    
    if any(strcmp(workspace_list, wksp))
        disp(strcat(' - ', wksp))
    else
        error(strcat('Fatal error: Workspace "', wksp, '" specified in chosen_workspaces not found in ', workspace_directory))
    end
end

%% Save any variables that may be overridden by workspace loading
remove_INH_control_temp = remove_INH_control;
remove_extra_controls_temp = remove_extra_controls;
remove_bad_treated_temp = remove_bad_treated;

%% Concatenate and set up workspaces

%~~~~~if more than one workspace was selected~~~~~%
if length(chosen_workspaces) > 1
    
    % turn into string for later
    chosen_workspace = "";
    for i = chosen_workspaces
        wksp = i{1};
        chosen_workspace = strcat(chosen_workspace, extractBefore(wksp,"."), " ");
    end
 
    % external function to combine workspaces
    final_data_table = combine_workspaces(chosen_workspaces,workspace_directory);
    
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
    
    if ~any(strcmp(final_data_table.Properties.VariableNames, 'DRUG_EXP'))
        % create drug-experiment column
        drug_exp = cellstr(strcat(final_data_table.DRUG,"_",final_data_table.EXP));
        final_data_table.DRUG_EXP = drug_exp;

        % reorganize table with new column
        final_data_table = [non_numeric final_data_table(:,end) numeric_data];
    end
    
    % Renormalize/zero_mean and calculate metadata variables for new
    % concatenated table
    merged_table_setup
    
    % Indicate that we've included multiple workspaces:
    conditions = true;
    
%~~~~~if only one workspace was selected~~~~~%
else
    % get name without .mat extension for later
    chosen_workspace = extractBefore(chosen_workspaces{1},".");
    % create full filepath
    fpath = strcat(workspace_directory,'/',chosen_workspaces{1});
    % load chosen workspace
    load(fpath)
    
    % there is only one condition/workspace, not multiple
    conditions = false;
end

%% Rewrite saved variables
remove_INH_control = remove_INH_control_temp;
remove_extra_controls = remove_extra_controls_temp;
remove_bad_treated = remove_bad_treated_temp;
clear remove_INH_control_temp remove_extra_controls_temp remove_bad_treated_temp

%% If TVN already perfomed, use pca_whitened_table ~~What if TVN has been merged with non-TVN?
%% If TVN already perfomed, use pca_whitened_table ** might need to sort pca_whitened too **

if ~exist('after_TVN','var')
    after_TVN = false; 
end
if after_TVN
    % set final_data_table to pca_whitened_table 
    final_data_table = pca_whitened_table;
end

%% Remove INH_control and other controls / bad treateds if selected
if remove_INH_control
    disp('Removing INH controls from final_data_table')
    final_data_table(contains(final_data_table.DRUG, 'INH_control'), :) = []; %remove the rows
    
end

if remove_extra_controls
    disp('Removing extra controls from final_data_table')
    controls = {'EtOH', 'MeOH', 'water', 'NaOH'};
    final_data_table(ismember(final_data_table.DRUG, controls), :) = []; %remove the rows
    clear controls
end

% TODO: this should account for IMI_dose etc. (use contains())
% Also: bad_treated should be in a metadata workspace probably
if remove_bad_treated 
    disp('Removing IMI, Nal, Nig, Dau')
    bad_treated = {'IMI', 'Nal', 'Nig', 'Dau'};
    final_data_table(ismember(final_data_table.DRUG, bad_treated), :) = []; %remove the rows
    clear bad_treated
end 

%%% Reinit after removal
if remove_INH_control || remove_extra_controls || remove_bad_treated % only do this once if either is triggered
    disp('Reinitializing metadata after drug removal')
    merged_table_setup %recreate metadata
end
%% Sort final_data_table by MoA
% sort by MoA so that drugs with similiar MoA are listed near each other
final_data_table = sort_by_column_name(final_data_table, 'ID');

%% Get a list of drugs and experiments

% get list of drugs from table
drugs = unique(final_data_table.DRUG,'stable')';
try
    experiments = unique(final_data_table.EXP, 'stable');
catch
    disp("no EXP column")
end
%% Initialize choice_variable and choice_extension
% If multiple workspaces, drugs will be denoted by drug_exp
% If only one workspace, drugs will be denoted by drug
%
% choice_variable -> d
%
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