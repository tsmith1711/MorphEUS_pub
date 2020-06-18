%% knn helper script 
% designed to be used by a bigger knn script. cuts down on repetitive code
% need to set up knn_data and knn_drugs ahead of time

% could maybe make a function? Inputs
% knn_data
% knn_drugs 
% d_metric - default to cosine 
% large_group - default to true
% create_graph - default to false
% create_cm -default to off 

% outputs would be
% pct_correct
% cm object if cm was created

if exist('cos_distance','var')
    if cos_distance
        d_metric = 'cosine';

    else
        d_metric = 'Euclidean';
    end
end 
%% create some defaults if not declared earlier 
% where to save figures 
if ~exist('figure_folder_path','var')
    figure_folder_path = "./figures/";
    disp('defaulted figure_folder_path to ./figures/');
end 

% if we're looping through and doing bayesian stuff
if ~exist('save_bayes_results','var')
    save_bayes_results = false;
end 

% distance metric to usee
if ~exist('d_metric','var')
    d_metric = 'cosine';
    disp('defaulted d_metric to cosine')
end

% if efflux should be in the lipid category
if ~exist('efflux_in_lipid','var')
    efflux_in_lipid = true;
    disp('defaulted efflux_in_lipid to true')
end

% if running a control confusion 
if ~exist('confuse_controls','var')
    confuse_controls = false;
end

% boolean whether to create graph
if ~exist('create_graph','var')
    create_graph = false;
    disp('defaulted create_graph to false')
end

 % set to on to create, off to not create (will be invisible with off)
if ~exist('create_cm','var')
   create_cm = 'off'; 
   disp('defaulted create_cm to off')
end

% which settings to use -- more general settings or more specific settings
    % true: use only the large 5 groups
% false: break into finer groups 
if ~exist('large_group','var')
    large_group = true;
    disp('defaulted large_group to true')
end 

% boolean: whether to do second nearest neighbors 
if ~exist('nn2','var')
   nn2 = false; 
   disp('defaulted nn2 to false')
end

% if knn was performed after pca, assume false if it doesn't exist 
if ~exist('knn_after_pca','var')
        knn_after_pca = "";
end 

%% knn color setup

knn_color_setup

%% do knn calculations 
% will rewrite this with the neighbors 
nearest_neighbors = knn_drugs;

% preallocate space to store distances
knn_distances = zeros(length(knn_drugs),1);

num_neighbors = 1;
% if doing second nearest neighbors 
if nn2
    nearest_neighbors2 = knn_drugs;
    knn2_distances = zeros(length(knn_drugs),1);
    num_neighbors = 2;
end 

for i = 1:size(knn_data,1)
    % get individual drug
    individual_row = knn_data(i,:);

    % get remainder
    remaining_rows = knn_data;
    remaining_rows(i,:) = [];

    % get the drugs for 'remaining_rows'
    remaining_drugs = knn_drugs;
    remaining_drugs(i,:) = [];

    % do knn
    [idx,distances] = knnsearch(remaining_rows,individual_row,'Distance',d_metric,'K',num_neighbors);
    knn_distances(i) = distances(1); 
    nearest_neighbors(i) = remaining_drugs(idx(1));
    if nn2
         knn2_distances(i) = distances(2);
         nearest_neighbors2(i) = remaining_drugs(idx(2));
    end 
end 

% store data in results table 
if nn2
results_table = table(knn_drugs,nearest_neighbors,knn_distances,nearest_neighbors2,knn2_distances, ...
        'VariableNames',{'DRUG','NEIGHBOR','DISTANCE','NEIGHBOR2','DISTANCE2'});
else
    results_table = table(knn_drugs,nearest_neighbors,knn_distances,'VariableNames',{'DRUG','NEIGHBOR','DISTANCE'});
end


%% ADD THING HERE FOR BAYSIAN STUFF 

drug_list = results_table.DRUG;
neighbor_list = results_table.NEIGHBOR;
if nn2
    neighbor2_list =results_table.NEIGHBOR2;
    
end 

% if doing baysean stuff, need to store results in a separate var

if save_bayes_results
    for i = 1:length(drug_list)
        drug = drug_list{i};
        neighbor = neighbor_list{i};
        bayes_struct.(drug).(neighbor) = bayes_struct.(drug).(neighbor)+1;
        if nn2
            % if we are alsol looking at second nearest neighbors
            neighbor2 = neighbor2_list{i};
            bayes_struct.(drug).(neighbor2) = bayes_struct.(drug).(neighbor2)+1;
            
        end 
    end 
    disp("saved bayes results")
end 

%% Create confusion matrix -- setup
% get moas for each drug -> done in external function now!! wow! 

drug_moas = change_to_moa(drug_list,all_categories);
%drug_moas = drug_list;
neighbor_moas = change_to_moa(neighbor_list,all_categories);
%neighbor_moas = neighbor_list;

% for control confusion
if confuse_controls
    drug_moas = drug_list;
    neighbor_moas = neighbor_list; 
end 

%% Create confusion matrix -- graph  
% confusion matrix figure will only show if user asked, otherwise it will
% be invisible (spooky)
f = figure('visible',create_cm);
% generate confusion matrix with row summaries 
cm = confusionchart(drug_moas,neighbor_moas, 'RowSummary','row-normalized');

% create title for figure
cm_tit = strcat("confusion matrix -- ", d_metric, " ", num2str(numvar), " variables ",tit_add);
title(cm_tit)

% use external function to get percent correct
pct_correct = cm_error(cm);

% display success percent in terminal window 
if strcmp(create_cm,"on")
    disp(strcat(num2str(pct_correct), "% for ", num2str(numvar), " variables"))
    % save figure with figure information
    ver_figsave(f,strcat(figure_folder_path,cm_tit),starting_vars,unique(knn_drugs),...
        strcat(num2str(pct_correct), "% for ", num2str(numvar), " variables"));
else
    % close the figure if it was invisible 
    close(f) 
end
%% Create graph

if create_graph 
    if nn2
        knn_digraph =  digraph([results_table.DRUG; results_table.DRUG],...
            [results_table.NEIGHBOR; results_table.NEIGHBOR2], ...
            [results_table.DISTANCE; results_table.DISTANCE2]);
    else
        knn_digraph =  digraph(results_table.DRUG,results_table.NEIGHBOR,results_table.DISTANCE);
    end 
    
    % get colors for each node 
    node_colors = assign_knn_color(knn_digraph,color_struct,all_categories);

    % create figure
    knn_fig = figure;
    knn_graph = plot(knn_digraph,'Layout','force','LineWidth',2.5,'NodeColor',node_colors);
    knn_axes = gca;
    
    
    % color based on weight
    knn_graph.EdgeCData= knn_digraph.Edges.Weight;

    % change background color
    set(knn_axes,'Color',[.8 .8 .8])
    
    % create title based on information 
    if strcmp(knn_type,'median')
        if knn_with_apply
            tit = strcat(knn_after_pca," - medians ", addedDrugStr," applied ",d_metric, " ", num2str(numvar)," variables ");
        else  
            % if we made this one by randomly swapping labels, want to show
            % that in the title
            try
                tit = strcat(knn_after_pca," - medians ",d_metric, " ", num2str(numvar)," variables ",tit_add);
            catch
                tit = strcat(knn_after_pca," - medians ",d_metric, " ", num2str(numvar)," variables ");
            end 
        end
    else
        if knn_with_apply
            tit = strcat("KNN after PCA ", addedDrugStr," applied ",d_metric, " ",num2str(numvar)," variables ");
        else
            tit = strcat("KNN after PCA ",d_metric);
        end
    end
    
    
    % add title 
    title(knn_axes, tit,'Interpreter','none')
    
    % add colorbar
    knn_cb = colorbar(knn_axes);
    knn_cm = colormap(knn_fig, map);
    
    % resize image
    set(knn_fig, 'Position', [100, 100, 800, 700])
    
    % save figure with external function 
    ver_figsave(knn_fig,strcat(figure_folder_path,tit),starting_vars,unique(knn_drugs),...
        strcat(num2str(pct_correct), "% for ", num2str(numvar), " variables"))
end 