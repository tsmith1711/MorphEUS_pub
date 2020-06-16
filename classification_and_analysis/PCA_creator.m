%% PCA Analysis with user prompts for workspace, number of variable, and drug selection
% DESCRIPTION:
%     Prompts users to obtain the following variables
%         - numvar
%         - chosen_workspaces
%         - chosen_drugs
%         - drugs_to_apply
% 
%     Having declared these variables, runs the following scripts in this
%     order:
%         - load_workspaces
%         - PCA_analysis
%
% PREREQUISITES: none
% 
% OUTPUT: ~fill this out~

%% Clear workspace before start
%~~ Comment out to disable this feature ~~%
%clear

%% Add Paths
add_all_paths 


%% Prompt for figure workspace

%figure_folder_path = inputdlg("Folder to save figures in:", "Figure Path", [1 35], "./figures/");

figure_folder_path = "./figures/";
%% Prompt for PCA_analysis settings

% Scripts like get_existing_feature_set set these variables, and we want
% them to be respected in the gui selection. If they're no set there
% though, they still need a default value.
if ~exist('do_feature_reduction','var')
    do_feature_reduction = true;
end

if ~exist('remove_extra_controls','var')
    remove_extra_controls = true;
end

disp('bad treated = IMI, Dau, Nig, Nal');
variables_to_create = {'disc', 'do_feature_reduction', 'include_DMSO', 'show_clustergram', 'show_pca', 'plot_3', 'plot_3_DMSO', 'plot_before_tvn', 'remove_INH_control', 'remove_extra_controls', 'remove_bad_treated', 'prompt_for_one_dose', 'use_medians_only'};
default_values =      [true,   do_feature_reduction,      false,            false,            true,       false,    false,       false,            true,                 remove_extra_controls,     true,               false,                  false];

target_values = checkboxList("Select PCA_analysis Settings", variables_to_create, default_values);

initialize_variables %Initialize variables_to_create as target_values

if use_medians_only
    numvar = 25;
    
end 

%% make this global to make it accessible inside ver_figsave()
global chosen_workspaces

%% Select number of variables to use

if do_feature_reduction
    prompt = {'Enter number of variables:'};
    dlgtitle = 'Variable Number';
    dims = [1 35];
    definput = {'94'};
    numvar = inputdlg(prompt,dlgtitle,dims,definput);
    numvar = str2num(numvar{1});
else
    if exist('pca_whitened_table','var')
        numvar =  sum(varfun(@isnumeric,pca_whitened_table,'OutputFormat', 'uniform'));
    else
        
    end 
end 
%% Find list of workspaces, prompt user for choice, and load chosen
if ~exist('workspace_directory', 'var')
    workspace_directory = './workspaces';
    disp("Defaulted workspace_directory to './workspaces' in load_workspaces")
end

workspace_list = get_workspaces(workspace_directory);
if isempty(workspace_list)
    error(strcat("No workspaces found in ", workspace_directory))
end
% If we are operating from an existing feature set, automatically select
% the workspaces that were used to calculate that set
if exist('overall_vars', 'var')
    try
        default_selection = find(ismember(workspace_list, overall_vars.WORKSPACE));
    catch
        disp("Note: overall_vars does not have a WORKSPACE field")
        default_selection = [];
    end 
else
    default_selection = [];
end

% prompt the user to select 
[indexes_chosen, any_chosen_tf] = listdlg('ListString',workspace_list,...
    'Name',"Select workspace",'ListSize',[380 380],...
    'InitialValue', default_selection);

% if the user did not select anything, exit (there is a catch for this in load_workspaces too)
if ~any_chosen_tf
    return
end

% use the index to get the list of user-selected workspaces
chosen_workspaces = workspace_list(indexes_chosen);

%% run the script to load the workspaces
load_workspaces

%% From the list of drugs, prompt user to select which drugs to use and apply

% If we are operating from an existing feature set, automatically select
% the variables that were used to calculate that set
if exist('overall_vars', 'var')
    default_selection = find(ismember(choice_variable, overall_vars.DRUGS));
else
    default_selection = [];
end

suggestion = {'Mer','Amp','Ctax','INH','EMB','ETA','IMI','Van','Cyc','Del',...
'Lev','Mox','Clz','MIT','Olf','Kan','Amk','Cam','Cla','Dox','Gent',...
'Strep','Tet','Lin','Pre','CCCP','Cer','Mon','Nig','Thi','RifT','BDQ',... 
'RIF','THL','water','Untreated'};

if exist('suggestion','var')
    default_selection = find(ismember(choice_variable,suggestion));
    
end 

% For extra_resolution_confusion
if prompt_for_one_dose
    select_one_dose
end

% Prompt user to select a list of drugs
[indexes_chosen, any_chosen_tf] = listdlg('ListString',choice_variable,...
    'Name',"Select Drugs",'ListSize',[220 380],'PromptString',"Select Drugs",...
    'InitialValue', default_selection); 

% if the user did not select anything, exit
if ~any_chosen_tf
    return
end

% Get a list of the drugs chosen
chosen_drugs = choice_variable(indexes_chosen);

%prompt the user to select drugs to apply if they want to 
[indexes_chosen, any_chosen_tf] = listdlg('ListString',choice_variable,'Name',"Select Drugs to apply",'ListSize',[220 380],...
    'CancelString','No Selection'); 

% Get a list of drugs to apply

drugs_to_apply = choice_variable(indexes_chosen);

%% Run PCA Analysis

PCA_analysis