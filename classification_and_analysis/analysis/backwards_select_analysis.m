%% backwards selecvt analysis
% all of the parts of backwards select without the prompts to load or save

% set this as as safety 
if ~exist('apply','var')
    apply = false;
end 

if ~exist('no_TVN','var')
    no_TVN = false; 
end 

%% Start backwards select 
% start off by getting original percet 
disp(strcat(num2str(pct_correct), "% for ", num2str(numvar), " variables"))

% this is what we are looking to beat
baseline_score = pct_correct;

% minumum variable number to run until 
minvar = 22; 

% set keep_going to true to get while loop started
keep_going = true;

% save all of the best variable iterations
overall_vars = struct; 

% use external function to get git info
git_info = getGitInfo();

% can save hash info in overall vars
try
    overall_vars.hash =  git_info.hash; 
catch
    % have a catch in case repository wasn't cloned 
    disp("unable to to save git hash information")
end 

% store the drugs used to create this in overall_vars 
overall_vars.DRUGS = unique(final_data_table.DRUG)';
% hard coding this for now b/c we're having problems 
% remove_after_TVN = true;
if remove_after_TVN
    % not actually using untreated, want to remove from .DRUG for clarity
    all_drugs_from_final = unique(final_data_table.DRUG)';
    unt_loca = find(strcmp(all_drugs_from_final,"Untreated"));
    all_drugs_from_final(unt_loca) = [];
    
    overall_vars.DRUGS = all_drugs_from_final;
end 
% save the workspace used
overall_vars.WORKSPACE = chosen_workspaces;

disp(strcat("distance metric is ", d_metric))

if large_group
    disp("using broader groups")
else
    disp("using finer groups")
end


% use randomize_labels function on the drug column in final_data_table
if with_randomized_labels
    
    % determine the randomized order
    rand_order = randperm(length(final_data_table.DRUG));

    disp("random!")
    overall_vars.RANDOM = true;
end 

% store in the overall_vars if joint profile or not 
if do_joint_profile
    overall_vars.JOINT_PROFILE = true;
else
    overall_vars.JOINT_PROFILE = false;
end 

% once the max value is not going up anymore, we stop this loop
disp("Beginning backwards select....")
disp("value of remove_after_TVN")
disp(remove_after_TVN)

while keep_going
    
    % want to make an empty % error thing to store
    all_pcts = zeros(length(starting_vars),1);
    
    % reset numvar 
    numvar = length(kept_vars);
    
    % get title for structure
    field_title = strcat("vars_",num2str(numvar),"_pct_",num2str(floor(baseline_score)));
    
    % save current number of variables in structure 
    overall_vars.(field_title) = kept_vars;
    
    %%% might also want to create a confusion matrix and save its handle in
    %%% overall_vars? 
    
    %% Now want to loop through, removing one variable at a time from starting_vars and calculating error 
    for p = 1:length(kept_vars)
        % reset to the full list of variables
        starting_vars = kept_vars;
        % get the name of the var removed
        var_removed = starting_vars{p};

        % remove a variable from the list
        starting_vars(p) = [];
        
        % perform TVN with newly selected features
        if no_TVN
            pca_whitened_table = final_data_table; 
        else
            TVN_transform
        end 
        % if applied drug, remove here
        if apply
            
            for qi = addedDrug
                added_drug_1 = qi{1};
                applied_inds = find(strcmp(pca_whitened_table.(choice_extension),added_drug_1));
            % remove!!!
                pca_whitened_table(applied_inds,:) = [];
            end 
 
        end 
        
        % if we want to get rid of the untreated we gotta do it every time
        %% Remove after TVN
        % forcing this true for now but change later
        if remove_after_TVN
            drugs = unique(pca_whitened_table.DRUG)';

            drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

            for i = drugs
                drug = i{1};
                drug_indexes.(drug) = find(strcmp(pca_whitened_table.DRUG, drug)).';
            end

            pca_whitened_table(drug_indexes.Untreated,:) = [];

            drugs = unique(pca_whitened_table.DRUG)';

            drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

            for i = drugs
                drug = i{1};
                drug_indexes.(drug) = find(strcmp(pca_whitened_table.DRUG, drug)).';
            end

        end 
        %% with randomized labels 
        if with_randomized_labels
            % go through each drug and reorder in the random order defined 
            randomized_drugs = pca_whitened_table.DRUG; 
            for w = 1:length(randomized_drugs)
                randomized_drugs(w) = pca_whitened_table.DRUG(rand_order(w));
            end 

            pca_whitened_table.DRUG = randomized_drugs;

        end 
        %% now time for the analysis 
        %find numeric columns of final table
        numeric_final_data_cols = varfun(@isnumeric,pca_whitened_table,'OutputFormat', 'uniform');

        %find column names
        numeric_final_col_names = pca_whitened_table.Properties.VariableNames(numeric_final_data_cols);
        
        if plot_before_tvn
             %find numeric columns of final table
            numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');
            
            numeric_final_col_names = final_data_table.Properties.VariableNames(numeric_final_data_cols);
            
            zero_mean_table = normalize_and_zero_mean(final_data_table,numeric_final_data_cols);

            table_to_zero_mean = final_data_table(:,starting_vars);

            table_to_zero_mean = horzcat(non_numeric,table_to_zero_mean);
            
            numeric_final_data_cols = varfun(@isnumeric,table_to_zero_mean,'OutputFormat', 'uniform');

            %find column names
            numeric_final_col_names = table_to_zero_mean.Properties.VariableNames(numeric_final_data_cols);

            zero_mean_table = normalize_and_zero_mean(table_to_zero_mean,numeric_final_data_cols);
            
            pca_whitened_table = zero_mean_table;
        end 
        % get numeric data
        ndata = pca_whitened_table{:,numeric_final_data_cols};
        
        % non numeric data
        non_numeric = pca_whitened_table(:,~numeric_final_data_cols);
        
        %% perform PCA calculation 
        [coeff,scores,pcvars] = pca(ndata);

        % put data in a table for better organization
        scores_table = horzcat(non_numeric,array2table(scores)); 

        drugs = unique(scores_table.DRUG)';

        % confirm drug_indicies 
        drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

        for i = drugs
            drug = i{1};
            drug_indexes.(drug) = find(strcmp(scores_table.DRUG, drug)).';
        end

        %% perform knn 
        if after_pca
            % if after pca, want to use scores table
            [drug_medians,drugs] = table_median(scores_table,choice_extension);
            overall_vars.AFTER_PCA = true;
        else
            % if not, want to use pca_whitened_table
            [drug_medians,drugs] = table_median(pca_whitened_table,choice_extension);
            overall_vars.AFTER_PCA = false;
        end
        

        % use medians as knn_data
        knn_data = drug_medians;
        knn_drugs = drugs';
        
        % if we're doing control confusion, needs to be like make_full
        if confuse_controls
            knn_data = scores;
            knn_drugs = scores_table.(choice_extension);    
        end 

        % run helper script to create confusion matrix (make sure things are
        % set to not show or else you will have A Time)
        knn_helper

        all_pcts(p) = pct_correct; 
        % now we have pct_error I guess just disp it for now yeeeet
     %  disp(strcat(num2str(pct_correct), "% for ",var_removed," removed"))
    end 
    
    % now want to find the max value in all_pcts
    [max_value] = max(all_pcts);
    
    % want to make sure we get all of the possible max indicies
    max_idx = find(all_pcts == max_value);

    %  if there is more than one location that has this max percent, want to
    %  only take one (and take the one that comes last)
    if length(max_idx) > 1
        max_idx = max_idx(end); 
    end 
    
    % if the max_value of the % succcess is better than the baseline score,
    % we like that and want to go with that new variable set
    % going to change this to 20 vars because generally don't trust
    % anything that low anyway, might as well speed things up
    % max_value >= baseline_score &&
    if length(starting_vars) > minvar 
        
      % want to throw this out from starting_vars and here we go ~again~
        kept_vars(max_idx) = []; 
        baseline_score = max_value;
        disp(strcat(num2str(baseline_score), "% and ", ...
            num2str(length(kept_vars)), " variables"))
        
    else
        keep_going = false;

        numvar = length(kept_vars);
        
    end 
end 
