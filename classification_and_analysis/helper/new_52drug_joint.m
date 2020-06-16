%% Script to create 025x-3x 52 drug joint profile
% Just going to get a rough working process of making a joint profile with
% the latest workspace so Michaela can work with this later

%% Discussion
% There are two "joint profile" pipelines:
%  - the "easyDoseResponse/doseResponse_analysis" pipeline, which:
%       - Extracts dose from the underscores after drugs, and then merges
%         drugs across a selcected amount of doses
%  - the "easyConditions/conditions_analysis" pipeline, which:
%       - Extracts "conditions" from the EXP column of the table, and joins
%         drugs across the selected EXPs
%
% These two pipelines were needed to deal with the original doseResponse
% and timecourse data, but in this instance, we could use either.
%
% DECISION: We're going to take the dose route, as I think it will make
% this a bit more consistent/robust. Either would work, though— and on
% principal I think working with an extra col rather than underscores is a
% bit better of an approach

%% Settings for analysis

if ~exist('bayesian','var')
   bayesian = false; 
end

% apply on 24hour timepoint workspace
apply_timepoint_drugs = false; %ALSO MUST COMMENT OUT CLUSTERGRAM SECTION OF PCA_analysis (idk lol it was throwing an error)

do_joint_profile = true;

% use defualt if did not use get exisitng feature set 
if ~exist('use_default','var')
    use_default = true;
end 

%PCA
show_pca = true;
plot_3 = false;
plot_3_DMSO = false;
show_clustergram = false;

%KNN
after_pca = true;
create_cm = 'on';
create_graph = true;
make_median = true;
make_full = false;
random_compare = false;

% set the rng seed for reproducibility
rng_val = 2 ;

disp("rng seed set to " + num2str(rng_val))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature Selection Options %

% default is to not do feature reduction and to use exisitng set

if ~exist('do_feature_reduction','var')
    
    do_feature_reduction = false; 
end 

if do_feature_reduction
    
    % %% MRMR ONLY -- UNCOMMENT TO RUN BACKWARDS SELECT -- (also uncomment the
   % clear at the beginning)
    numvar = 101;
    % for the bayes stuff, can only do 79 vars b/c 80 untreated, but will
    % start small for tests
    if bayesian
        if debug_mode 
            numvar = 31;
        else
            numvar = 79;
        end 
    end
    do_feature_reduction = true;

else 
    if use_default
%% USE EXISTING FEATURE SET 
        % bypasses get_existing_feature_set 
        overall_variables_path = './feature_sets/overall_vars_noUnknown_doseResponse_8-16-19.mat';
        variable_set_name = 'vars_77_pct_82';
        get_existing_feature_set_no_gui
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end 



%% Load the main workspace

% Drug removal
remove_INH_control = true; %This has to be true bc the underscore in the name messes joint profiles up

% load the workspace
% changed to use old DMSO for now
% chosen_workspaces = {'full_025x_and_3x_workspace_old_DMSO.mat'};
%chosen_workspaces = {'full_025x_and_3x_workspace.mat'};
chosen_workspaces = {joint_wksp};
load_workspaces

%% Split the tables so that we can manipulate them separately (Mwahaha)

% Find the 025x table
rows_025x = strcmpi(final_data_table.EXP, 'x025'); % Some of the EXPS are 'x025', and others 'X025'
final_data_table_025x = final_data_table(rows_025x, :);

% Find the 3x table
rows_3x = strcmp(final_data_table.EXP, 'x3');
final_data_table_3x = final_data_table(rows_3x, :);

% Find the 3x untreateds
rows_025x_unt = strcmp(final_data_table_025x.ID, 'Untreated');
untreateds_025x = final_data_table_025x(rows_025x_unt, :);

% Make sure we didnt miss any
cum_size = (size(final_data_table_025x, 1) + size(final_data_table_3x, 1));
if cum_size ~= size(final_data_table, 1)
    error("You missed some")
end


%% Cool table fixes
% Eventually, you'll probably want to do all these and save them in the
% workspace

% Set DRUG to drug_dose
final_data_table_3x.DRUG = strcat(final_data_table_3x.ID, '_3x');

%025x fixes
% Set drug to DRUG_dose
% This is already done, but this accoutns for some of the errors (cerold)
% We also want our controls to have the _025x now, as that will join them
% laterally with the 3x  (untreated -> untreated_025x)
final_data_table_025x.DRUG = strcat(final_data_table_025x.ID, '_025x');


%% Find the ids of each table, and the ones they have in common

ids_3x = unique(final_data_table_3x.ID).';
ids_025x = unique(final_data_table_025x.ID).';

% Find the common ids
common_ids = intersect(ids_3x, ids_025x);

removed_025x_ids = ids_025x(~ismember(ids_025x, common_ids));
removed_3x_ids = ids_3x(~ismember(ids_3x, common_ids));

% Let's talk about our feelings
if ~isempty(removed_025x_ids)
    disp("Removed the following ids from the 025x workspace:")
    showcell(strcat(" - ", removed_025x_ids));
end

if ~isempty(removed_3x_ids)
    disp("Removed the following ids from the 3x workspace:")
    showcell(strcat(" - ", removed_3x_ids));
end


%% Join our hands together in love and matrimony
final_data_table = [final_data_table_025x; final_data_table_3x];

% Set the EXP to the same thing so the joint profile allows them to canoodle
final_data_table.EXP(:) = {'drug52_doseResponse'};
%% Strike fear into the eyes of the enemy and only allow those we chose to survive

%Extract only the ids that the two doses have in common
final_data_table = final_data_table(ismember(final_data_table.ID, common_ids), :);


