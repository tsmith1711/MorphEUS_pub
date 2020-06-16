%% Script to perform KNN on data outputted by PCA_analysis
% This script uses existing variables in the workspace to perform KNN
% and output grahps
%
% PREREQUISITES:
% - ** Must have run PCA_analysis.m and all of its prerequisites **
%       - Most important variable is scores_table
% - All required settings in section below
% - Any desired optional variables listed below
%
% OUTPUT:
% - 
% 
% DEPENDENCIES:
% - knn_helper ***
%
%% Required setting variables
%%% Must include all of these variables in the workspace before running %%%
%
% - make_median - boolean whether to make a plot of replicate medians
% - make_full   - boolean whether to make a plot with each replicate as a point
% d_metric, create_graph, create_cm, large_group

%% Optional settings
%%% Optionally include these variables in the workspace before running if %
%%% desired. If not present, they will default to the specified value     %
%
% - random_compare - (default: false) - make a randomized label plot to compare
% - tit_add        - (default: "")    - 
%% Add paths
addpath('./functions')
addpath('./helper/')

%% Generate defaults for optional settings
if ~exist('random_compare', 'var')
    random_compare = false;
    disp('defaulted random_compare to false');
end

if ~exist('tit_add', 'var')
    tit_add = "";
    disp('defaulted tit_add to ""');
end

if ~exist('confuse_controls','var')
   confuse_controls = false;  
end

if confuse_controls
    make_median = false;
    make_full = true;
    
end 

disp("begining of knn_analysis, value of remove_after_TVN")
disp(remove_after_TVN)

%% Validate required settings

%From PCA_analysis -- just check for these two; they're the most important
if ~exist('scores_table', 'var')
    error('Fatal error: Must run PCA_analysis before KNN_analysis; missing scores_table')
end

if ~exist('choice_variable', 'var')
    error('Fatal error: Must run PCA_analysis before KNN_analysis; missing choice_variable')
end

%Extra required variables for this script
if ~exist('make_median', 'var')
    error('Fatal error: Missing required settings variable make_median in KNN_analysis')
end

if ~exist('make_full', 'var')
    error('Fatal error: Missing required settings variable make_full in KNN_analysis')
end
  
%% -- KNN analysis begins here -- 

%% if doing full analysis with randomized labels (most likely just for backwads_select purposes)


%% Using each image after PCA
% using scores table--reassigning to new variable here to make it flexible
if make_full
    
    knn_type = 'full';
    
    knn_with_apply = false;
    
    knn_table = scores_table;
        
    knn_data = knn_table{:,numeric_final_data_cols};

    knn_drugs = knn_table.(choice_extension);

    knn_helper
end 
%% Using medians after PCA

if make_median 
    tit_add = "";

    knn_type = 'median';
    knn_with_apply = false;
    if after_pca
        % if after pca, want to use scores table
        knn_after_pca = "KNN after PCA";
        % use separate script to get median values for each set of drugs in
        % knn_table 
        [drug_medians,drugs] = table_median(scores_table,choice_extension);
    else
        % if before pca, want to use pca_whitened_table
        knn_after_pca = "KNN before PCA";
        % use separate script to get median values for each set of drugs in
        % knn_table 
        [drug_medians,drugs] = table_median(pca_whitened_table,choice_extension);
    end
  
    
    %~~~ Perform knn ~~~%
    % go through each drug, make a separate row for the indiviudal drug and
    % compare it against the remainder 

    % use medians as knn_data
    knn_data = drug_medians;
    knn_drugs = drugs';

    disp("list of drugs for knn in knn_analysis")
    disp(knn_drugs)
    % run helper script 
    knn_helper

    % compare to random labels 
    if random_compare
        disp("Rerunning with random labels")
        tit_add = " with random labels";
        % copy list of knn_drusg 
        random_drugs = randomize_list(knn_drugs);
        % assign knn_drugs to be the new randomly created set
        knn_drugs = random_drugs;
        knn_helper 
    end 
end


%% knn for applied -> individual

if apply
    if make_full
        knn_type = 'full';
        knn_with_apply = true;
        % using scores table--reassigning to new variable here to make it flexible
        knn_table = newspace_table;

        knn_data = knn_table{:,numeric_final_data_cols};

        knn_drugs = knn_table.DRUG;

        knn_helper
    end 
end 

%% create knn for applied 

if apply
    if make_median
        knn_type = 'median';
        knn_with_apply = true;

        % get median values with external function 
        [drug_medians,drugs] = table_median(newspace_table,choice_extension);
        
    %% Perform knn
        % go through each drug, make a separate row for the indiviudal drug and
        % compare it against the remainder 

        % use medians as knn_data
        knn_data = drug_medians;

        knn_drugs = drugs';

        knn_helper 
        
        disp("did applied knn")
    end  
end 