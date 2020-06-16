%% PCA Analysis Script
%DESCRIPTION:
%   functionality-only pca analysis for data outputted by bcp-microbej-processing
%   Adapted from the easyPCA script
%
% PREQUISITES:
% - Must have run load_workspaces.m
% - Must define all "input variables" detailed below in the workspace
%   before running
%
% OUTPUT:
% - PCA scores, transformation, and loadings (+ more) in workspace
% - PCA figures - either displayed/variables/files (?)
%
% DEPENDENCIES
%  - filter_drugs
%  - normalize_and_zero_mean()
%  - TVN_transform
%  - feature_reduce
%  - (list incomplete)
%
%% TODO
% - Incorporate / make modular backwards_select
% - unify option selection
% - unify output / display of figures and data

%% Required variables in existing workspace:
% - all output of merged_table_setup.m, most importantly: final_data_table, drugs,
%    numeric_cols (etc.)
% - chosen_drugs: Horizontal cell array of drug names to include in calculations
% - drugs_to_apply: Horizontal cell array of drug names to apply to PCA�do
%    not include, or leave empty to apply no drugs
% - numvar: number of variables to use

%% Optional settings
%%% Include these in the workspace before running
%%% These defaults are generated below

% --> disc (default: true) % discretize for feature reduction
% --> do_feature_reduce (default: true) % perform mRMR feature reduction

% whether or not to use DMSO in TVN (Ando Normalization) calculations
% set to true to include DMSO, false to just use Untreated
% MUST SELECT DMSO AS A DRUG IF THIS IS SET AS TRUE
% --> include_DMSO (default: false)

% --> show_clustergram (default: true) % booleaen, whether or not to show clustergram
% --> show_pca (default: true) % boolean, whether or not to show PCA
% --> plot_3 (default: true) % booleaen: whether to plot only 3 untreated points in pcas

%%% To add:
% graph saving + path config (interface with saving script)


%% Description of action * out of date *
% 0. numvar
% 1. Workspace
% 2. Drugs
% 3. Applied drug
% 4. Performs mrmr feature selection w/ number of features as chosen by
% user in the workspace selection/setup section 
% 5. Performs TVN on data using untreatd controls (and DMSO if user selected that option) 
% 6. Based on previous selections and user settings, will display PCAs and
% clustergrams 

%% Add all needed paths
addpath('./functions')  
addpath('./helper/')

% add paths for mrmr
addpath('./mRMR_0.9_compiled')
addpath('./mRMR_0.9_compiled/mi')

%% Generate defaults for missing user options

% discretize for feature reduction
if ~exist('disc', 'var')
    disc = true;
    disp('defaulted disc to true')
end

% whether to do feature reduction
if ~exist('do_feature_reduction', 'var')
    do_feature_reduction = true; 
    disp('defaulted do_feature_reduction to true');
end

% whether or not to use DMSO in TVN calculations
% set to true to include DMSO, false to just use Untreated
% MUST SELECT DMSO AS A DRUG IF THIS IS SET AS TRUE
if ~exist('include_DMSO', 'var')
    include_DMSO = false;
    disp('defaulted include_DMSO to false');
end

% booleaen, whether or not to show clustergram
if ~exist('show_clustergram', 'var')
    show_clustergram = false;
    disp('defaulted show_clustergram to false');
end

% boolean, whether or not to show PCA
if ~exist('show_pca', 'var')
    show_pca = false;
    disp('defaulted show_pca to false');
end

% boolean
if ~exist('no_TVN','var')
    no_TVN = false;
    
end 

% if not doing bayesian
if ~exist('bayesian','var')
   bayesian = false; 
end

% booleaen: whether to plot only 3 Untreated in pca
if ~exist('plot_3', 'var')
    plot_3 = false;
    disp('defaulted plot_3 to false');
end

% booleaen: whether to plot only 3 DMSO in pca
if ~exist('plot_3_DMSO', 'var')
    plot_3_DMSO = false;
    disp('defaulted plot_3_DMSO to false');
end

% boolean: whether to create a pca plot of the data before TVN
if ~exist('plot_before_tvn','var')
    plot_before_tvn = false;
    disp("defaulted plot_before_tvn to false")
end

% figure folder path
if ~exist('figure_folder_path','var')
    figure_folder_path = "./figures/";
    disp("defaulted figure_folder_path to ./figures")
end

if ~exist('remove_after_tvn','var')
    remove_after_tvn = false;
end 

%% Verify and initialize drug selections
filter_drugs

%% Re- normalize and zero-mean data after filtering drugs
% This helps get a better feature reduction

% if we've already done TVN, no need to zero mean
if after_TVN
    zero_mean_table = final_data_table; 
else 
    zero_mean_table = normalize_and_zero_mean(final_data_table, numeric_final_data_cols);
end

%% Feature reduce

if do_feature_reduction
    % perform feature reduction (separate script)
    feature_reduce %sets starting_vars to kept_vars
elseif exist('starting_vars','var')
    % use a predetermined user set 
    kept_vars = starting_vars;
    numvar = length(starting_vars);
else 
    %Just use all variables
    starting_vars = numeric_final_data_cols;
end

% this is somehow becoming false somewhere even if I want it to be true so
% look into that later I guess
%remove_after_TVN = true;

%% Perform TVN transform

% if there's an applied drug, add it back in
if apply 
    final_data_table = vertcat(final_data_table, drug_to_apply);
end 

% perform TVN transformation (if not done yet), else just implement
% starting variables chosen in previous step
if after_TVN
    % don't need to run this bad boy again if we've already done it
    % but still want to make sure we apply the starting/kept vars
    
    normalized_table = final_data_table;
    %find numeric columns of final table
    numeric_final_data_cols = varfun(@isnumeric,normalized_table,'OutputFormat', 'uniform');
    
    % get non numeric data
    non_numeric = normalized_table(:,~numeric_final_data_cols); 

    % get the variables that are kept
    kept_variables = normalized_table(:,starting_vars);

    normalized_table = horzcat(non_numeric,kept_variables);
    
    pca_whitened_table = normalized_table;
    
else
    if no_TVN
        pca_whitened_table = final_data_table;
    else 
        TVN_transform
        % if we forced user to select Untreated for the sake of the transform,
        % we should get rid of it here. 
        if remove_after_TVN

            drugs = unique(pca_whitened_table.DRUG)';

            drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

            for i = drugs
                drug = i{1};
                drug_indexes.(drug) = find(strcmp(pca_whitened_table.DRUG, drug)).';
            end

            pca_whitened_table(drug_indexes.Untreated,:) = [];
            disp("Removed Untreated after TVN transform")
            numeric_final_data_cols = varfun(@isnumeric,pca_whitened_table,'OutputFormat', 'uniform');
        end 
    end 
end 


%% Get numeric columns and check that indicies are correct

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,pca_whitened_table,'OutputFormat', 'uniform');

%find column names
numeric_final_col_names = pca_whitened_table.Properties.VariableNames(numeric_final_data_cols);

% get list of drugs in table
drugs2 = unique(pca_whitened_table.DRUG)';

drug_indexes = cell2struct(cell(1,length(drugs2)), drugs2, 2);

for i = drugs2
    drug = i{1};
    drug_indexes.(drug) = find(strcmp(pca_whitened_table.DRUG, drug)).';
end

% if we had applied a drug, we need to remove it from the table before
% doing PCA calculations
% also want to make a copy of the table 
if apply
    % choice_extension: DRUG if not combined workspace, DRUG_EXP if combind
    % workspace
    removal_var = unique(pca_whitened_table.(choice_extension))';
    
    remove_indexes = cell2struct(cell(1,length(removal_var)), removal_var, 2);

    for i = removal_var
        drug = i{1};
        remove_indexes.(drug) = find(strcmp(pca_whitened_table.(choice_extension), drug)).';
    end

    % make copy of table
    whitened_table_copy = pca_whitened_table;

    % now remove the added drug from the table, making 
    all_removed_inds = [];
    for i = addedDrug
        drug = i{1};
        all_removed_inds = [all_removed_inds remove_indexes.(drug)];
    end 
    
    pca_whitened_table(all_removed_inds,:) = [];
    
end 

%% Perform PCA calculations
ndata = pca_whitened_table{:,numeric_final_data_cols}; %numeric data from normalized table 

% non numeric data
non_numeric = pca_whitened_table(:,~numeric_final_data_cols);

% Columns of the scores variable corresponds to the principal components (from first to last, decreasing order)
% dataInPrincipalComponentSpace = ndata*coeff
[coeff,scores,pcvars] = pca(ndata);
vari = cumsum(pcvars)./sum(pcvars);
coeff=coeff';
pcvari = pcvars/sum(pcvars);

% put data in a table for better organization 
scores_table = horzcat(non_numeric,array2table(scores));

numeric_final_col_names_space = numeric_final_col_names;
%Alters the column names (convert underscores (_) to spaces)
for i=1:size(numeric_final_col_names,2)
    numeric_final_col_names_space(1,i) = strrep(numeric_final_col_names(1,i),'_',' ');
end

% find all principal components and percentages for axes labels 
pc_labels = cell(1,length(pcvari));
for i = 1:length(pcvari)
    pc_labels(i) = cellstr(strcat("PC:",num2str(round(i,2)), " ", (num2str(pcvari(round(i,2))*100)), "%"));
end

% on the chance that numeric cols changed b/c not enough data
scores_cols =  varfun(@isnumeric,scores_table,'OutputFormat', 'uniform');

%% PCA calculations before TVN 

if remove_after_TVN
    % remove untreated if user did not select them 
    % this is for before PCA data only !
    unt_rows = contains(zero_mean_table.DRUG,'Untreated');
    zero_mean_table(unt_rows,:) = [];
    % recalculate indicies

end 

before_tvn_data = zero_mean_table{:,numeric_final_data_cols}; %numeric data from normalized table 

% non numeric data
before_tvn_non_numeric = zero_mean_table(:,~numeric_final_data_cols);

% Columns of the scores variable corresponds to the principal components (from first to last, decreasing order)
[before_tvn_coeff,before_tvn_scores,before_tvn_pcvars] = pca(before_tvn_data);
before_tvn_vari = cumsum(before_tvn_pcvars)./sum(before_tvn_pcvars);
before_tvn_coeff=before_tvn_coeff';
before_tvn_pcvari = before_tvn_pcvars/sum(before_tvn_pcvars);

% put data in a table for better organization 
before_tvn_scores_table = horzcat(before_tvn_non_numeric,array2table(before_tvn_scores));

% find all principal components and percentages for axes labels 
before_tvn_pc_labels = cell(1,length(before_tvn_pcvari));
for i = 1:length(before_tvn_pcvari)
    before_tvn_pc_labels(i) = cellstr(strcat("PC:",num2str(round(i,2)), " ", (num2str(before_tvn_pcvari(round(i,2))*100)), "%"));
end

%% PCA Plot
drugs = unique(scores_table.DRUG)';

% confirm drug_indicies in case an applied drug was removed
drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

for i = drugs
    drug = i{1};
    drug_indexes.(drug) = find(strcmp(scores_table.DRUG, drug)).';
end

%% only plot 3 untreated
% confirm that untreated field does exist 
if plot_3 && any(contains(fieldnames(drug_indexes),'Untreated'))
    if conditions
        % if we have multiple conditions, neeed to get untreated from each
        % condition
        
        % get drug_experiments
        drug_exps = unique(scores_table.DRUG_EXP)';

        % get indicies for the drug-experiment column
        drug_exp_indexes = cell2struct(cell(1,length(drug_exps)), drug_exps, 2);

        for i = drug_exps
            drug = i{1};
            drug_exp_indexes.(drug) = find(strcmp(scores_table.DRUG_EXP, drug)).';
        end
        
        % find indicies of different types of untreated
        diff_untreated = contains(drug_exps,'Untreated');
        untreated_types = drug_exps(diff_untreated);
        
        % remove most untreeated, only plot 3
        for p = untreated_types
            drug = p{1};
            scores_table(drug_exp_indexes.(drug)(4:end),:) = [];
            
            % recalculate indicies 
            drug_exp_indexes = cell2struct(cell(1,length(drug_exps)), drug_exps, 2);

            for i = drug_exps
                drug = i{1};
                drug_exp_indexes.(drug) = find(strcmp(scores_table.DRUG_EXP, drug)).';
            end

        end
        
    else
        % remove most untreated, only plot 3
        scores_table(drug_indexes.Untreated(4:end),:) = [];
        before_tvn_scores_table(drug_indexes.Untreated(4:end),:) = [];
        
        % remove most untreated from pca_whitened_table as well (helps for the
        % later knn code)
        pca_whitened_table(drug_indexes.Untreated(4:end),:) = [];
        zero_mean_table(drug_indexes.Untreated(4:end),:) = [];
        
        % recalcualte indexes 
        drug_indexes = cell2struct(cell(1,length(chosen)), chosen, 2);

        for i = chosen
            drug = i{1};
            drug_indexes.(drug) = find(strcmp(scores_table.DRUG, drug)).';
        end
        
    end
end 

if plot_3_DMSO && any(contains(fieldnames(drug_indexes),'DMSO'))
   if conditions 
        % get drug_experiments
        drug_exps = unique(scores_table.DRUG_EXP)';

        % get indicies for the drug-experiment column
        drug_exp_indexes = cell2struct(cell(1,length(drug_exps)), drug_exps, 2);

        for i = drug_exps
            drug = i{1};
            drug_exp_indexes.(drug) = find(strcmp(scores_table.DRUG_EXP, drug)).';
        end
        
        % if DMSO exists, only want 3 of those as well 
        if isfield(drug_indexes,'DMSO')
            % find indicies of different types of untreated
            diff_DMSO = contains(drug_exps,'DMSO');
            DMSO_types = drug_exps(diff_DMSO);
            
            
             % remove most untreeated, only plot 3
            for p = DMSO_types
                drug = p{1};
                scores_table(drug_exp_indexes.(drug)(4:end),:) = [];

                % recalculate indicies 
                drug_exp_indexes = cell2struct(cell(1,length(drug_exps)), drug_exps, 2);

                for i = drug_exps
                    drug = i{1};
                    drug_exp_indexes.(drug) = find(strcmp(scores_table.DRUG_EXP, drug)).';
                end

            end
        end 
       
   else
       
        % if there are DMSO, only plot 3 as well
        if isfield(drug_indexes,'DMSO')
            % just nab a few DMSO to plot 
             scores_table(drug_indexes.DMSO(4:end),:) = [];
            % recalcualte indexes 
            drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

            for i = drugs
                drug = i{1};
                drug_indexes.(drug) = find(strcmp(scores_table.DRUG, drug)).';
            end
        end 
       
   end 
    
    
end

%% PCA figures

%creates distinct colors when plotted
c = distinguishable_colors(size(chosen,2),{'w','k'}); 

% get scores again -- will be same as before if no need to remove Untreated
scores = scores_table{:,scores_cols};

 %~~~~~~~~~~ 3D PCA Plot ~~~~~~~~~~~%

if show_pca
    
    % create title including var # and workspace used 
    tit = strcat(chosen_workspace,"-", num2str(sum(numeric_final_data_cols)), 'variables');
    
    % adjust title to indicate if DMSO was used in calculations
    if include_DMSO
        tit = strcat(tit, ' with DMSO');
    end
    
    % use external function to plot 
    % if user wanted a before tvn pca
    if plot_before_tvn
        before_tvn_hpca = pca_plot(before_tvn_scores_table,before_tvn_pc_labels,choice_extension);
        before_tvn_tit = strcat(chosen_workspace,"-", num2str(sum(numeric_final_data_cols)), 'variables before TVN');
        % set figure title
        title(before_tvn_tit,'Interpreter','none')
        % get save path title 
        figpath = strcat(figure_folder_path,before_tvn_tit);
        % save figure with metadata 
        ver_figsave(before_tvn_hpca,figpath,numeric_final_col_names,chosen)
    end 
    % after tvn pca (the usual one) 
    hpca = pca_plot(scores_table,pc_labels,choice_extension);
   

    % set figure title
    title(tit,'Interpreter','none')

    % get save path title 
    figpath = strcat(figure_folder_path,tit);
    % save figure with metadata 
    ver_figsave(hpca,figpath,numeric_final_col_names,chosen)

end 
%% Create Clustergram
% if user selected show clustergram
if show_clustergram
    
    % clustergram needs each row to have a unique value, therefore instead
    % of BDQ BDQ BDQ we need BDQ1 BDQ2 BDQ3
    yvalues = (scores_table.DRUG);
    for i = chosen
        drug = i{1};
        if conditions
            % use drug_exp indexes 
            drug_exps = unique(scores_table.DRUG_EXP);

            drug_exp_indexes = cell2struct(cell(1,length(drug_exps)), drug_exps, 2);

            for q = drug_exps
                drug = q{1};
                drug_exp_indexes.(drug) = find(strcmp(scores_table.DRUG_EXP, drug)).';
            end

            indicies = drug_exp_indexes;
        else
            indicies = drug_indexes;
        end
        for j = 1:length(indicies.(drug))
            k = indicies.(drug)(j);
            yvalues(k) = strcat(yvalues(k),num2str(j));
        end
    end
    % convert underscore to space for ease of reading
    yvalues = cellstr(strrep(yvalues,'_',' '));

    % set up colormap for clustergram 
    vec = [100; 87.5; 75; 62.5; 50; 37.5; 25; 12.5; 0];
    hex = ['#260000';'#610000';'#b00000';'#ff7575';'#feffff';'#87e2ff';'#00b0e9';'#006687';'#002b39'];
    raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
    N = 128;
    % get PC percents for labels
    colLabel1 = strcat("PC1 (",num2str(vari(1,1)*100),"%)");
    colLabel2 = strcat("PC2 (",num2str(vari(2,1)*100 - vari(1,1)*100),"%)");
    colLabel3 = strcat("PC3 (",num2str(vari(3,1)*100 - vari(2,1)*100),"%)");
    colLabels = cellstr([colLabel1 colLabel2 colLabel3]);
    map = interp1(vec,raw,linspace(100,0,N),'pchip');

    % make clustergram
    z = clustergram(scores(:,1:3),'RowLabels',yvalues,'ColumnLabels',colLabels,'Colormap',map);
    addTitle(z,tit)
    
    
    % want to do pre-TVN clustergram if user wanted that 
    if plot_before_tvn
        % before_tvn_scores_table,before_tvn_pc_labels

        colLabels = before_tvn_pc_labels(1:3);
        before_tvn_scores_cols = varfun(@isnumeric,before_tvn_scores_table,'OutputFormat', 'uniform');
        before_tvn_scores = before_tvn_scores_table{:,before_tvn_scores_cols};
        z = clustergram(before_tvn_scores(:,1:3),'RowLabels',yvalues,'ColumnLabels',colLabels,'Colormap',map);
        addTitle(z,"Before TVN")
    end 
end 

%% create PCA plot with applied drug -- only if user selected drugs to apply
%% Transform applied drugs into PCA space
if apply 
    
    normalized_table= whitened_table_copy;

    %find numeric columns of final table
    numeric_final_data_cols = varfun(@isnumeric,normalized_table,'OutputFormat', 'uniform');

    %find column names
    numeric_final_col_names = normalized_table.Properties.VariableNames(numeric_final_data_cols);

    ndata = normalized_table{:,numeric_final_data_cols}; %numeric data from normalized table 

    non_numeric = normalized_table(:,~numeric_final_data_cols);
    %apply eigevnectors
    newspace = ndata*coeff';
    
    % organize newspace into a table
    newspace_table = horzcat(non_numeric,array2table(newspace));

    drugStr = "";
    for i = chosen
        drug = i{1};
        drugStr = strcat(drugStr, drug, " ");
    end

    drugs = unique(normalized_table.DRUG)';

    %find row indexes for each drug
    drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

    for i = drugs
        drug = i{1};
        drug_indexes.(drug) = find(strcmp(normalized_table.DRUG, drug)).';
    end
end 

%% only plot 3 untreated and DMSO 
if apply
    if plot_3
        % just nab a few untreated to plot 
        newspace_table(drug_indexes.Untreated(4:end),:) = [];
        % recalcualte indexes 
        drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

        for i = drugs
            drug = i{1};
            drug_indexes.(drug) = find(strcmp(newspace_table.DRUG, drug)).';
        end
        

    end 
    if plot_3_DMSO
        
        % if there are DMSO, only plot 3 as well
        if isfield(drug_indexes,'DMSO')
            % just nab a few DMSO to plot 
             newspace_table(drug_indexes.DMSO(4:end),:) = [];
            % recalcualte indexes 
            drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

            for i = drugs
                drug = i{1};
                drug_indexes.(drug) = find(strcmp(newspace_table.DRUG, drug)).';
            end
        end 
        
    end 
end 

%% Create figure with applied drugs

% we don't need to go in here for bayesian stuff so let's just force us out
% and keep it easy breezy lemon squeezy
if apply && ~bayesian
    %creates distinct colors when plotted
    c = distinguishable_colors(size(drugs,2),{'w','k'}); 

    % get scores again -- will be same as before if no need to remove Untreated
    newspace = newspace_table{:,numeric_final_data_cols};
    
    addedDrugStr = "";
    for i = addedDrug
        drug = i{1};
        addedDrugStr = strcat(addedDrugStr, drug, " ");
    end
    
    if show_pca

        handpca2 = pca_plot(newspace_table,pc_labels,choice_extension);
        % set title
        tit = strcat(chosen_workspace,"- ", num2str(sum(numeric_final_data_cols))," variables with ", addedDrugStr," applied");
        title(tit,'Interpreter','none')
        
        % get figure path including the figure_folder
        figpath = strcat(figure_folder_path,tit);
        
        % save figure with metadata
        ver_figsave(handpca2,figpath,starting_vars,drugs,strcat(addedDrugStr," applied"))
            
    end 
    
    % create row titles for clustergram
    yvalues = (newspace_table.DRUG);
    for i = drugs2
        drug = i{1};
        for j = 1:length(drug_indexes.(drug))
            k = drug_indexes.(drug)(j);
            yvalues(k) = strcat(yvalues(k),num2str(j));
        end
    end

    % only create clusteregram if user wants it
    if show_clustergram
        % create clustergram
        z1 = clustergram(newspace(:,1:3),'RowLabels',yvalues,'ColumnLabels',colLabels,'Colormap',map);
    end
end 