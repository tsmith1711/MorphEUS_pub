%% to not repeat stuff, used by bayes_backwards_select
% content of the backwards select, 
% created in bayes_backwards_select:
% maxrun - max # of runs
% runct - current run number 
% do_joint_profile - whether or not this is a joint profile situation 
% suggestion - a list of suggested drugs to start with 


%%% need to add somewhere for random relabeling %%% 
%%% code from knn_helper on random analysis 
% copy list of knn_drusg 
% random_drugs = randomize_list(knn_drugs);
% % assign knn_drugs to be the new randomly created set
% knn_drugs = random_drugs;

% will begin by listing run number of total runs (maxrun) declared in
% bayes_backwards_select

disp(strcat("RUN NUMBER ", num2str(runct), " of ", num2str(maxrun)))

% want to start off by doing our initial mrmr to rank order, so need to
% start off with do_feature_reduction
do_feature_reduction = true;
% do not want to save bayes results after each run throughout
% backwards_select, so this must be off
save_bayes_results = false;

% numvar must be one fewer than the number of rows of untreated
% kept 
numvar = num_rows-1;
%% Set up workspace 

if do_joint_profile
    % same here for flattened and non flattened table -- at this point,
    % just puts the 3x and 025x in the same table
    new_34drug_joint
    
    %% prompt user selection of drugs
    % now that I've moved this around, could combine better with the
    % individual version--maybe do the picking later? 
    
    % if we're doing the split one, allow picking both high and low doses
    % so we want this to be drug and not id 
    if flattened_joint
        
        choice_variable = common_ids;
        choice_extension = 'ID';
        % if on later runs and applying condition data, need to make sure
        % to load those workspaces too 
        if runct > 1 && apply_timecourse_data
            load('./workspaces/all_timecourse_doseresponse.mat')
            
            table_combine_temp = struct;
            table_combine_temp.final_data_table = final_data_table;
            
            table_combine_temp.time = timecourse_data;
            
            big_data = column_validation(table_combine_temp);

            final_data_table = vertcat(big_data{:});
            
            choice_variable = unique(final_data_table.(choice_extension))'; 
                                
        end 
        
        
    else
        suggestion = [strcat(suggestion,"_025x"), strcat(suggestion,"_3x")];
        
        
        
        choice_variable = unique(final_data_table.DRUG)';
        choice_extension = 'DRUG';
        
    end 
    
    %% drug selection 
    % if chosen drugs already exists, skip this 
    if ~exist('chosen_drugs','var')
        if exist('overall_vars', 'var')
            default_selection = find(ismember(choice_variable, overall_vars.DRUGS));
        else
            default_selection = [];
        end

        % suggestion created in bayes_backwards_select 
        if exist('suggestion','var')
            default_selection = find(ismember(choice_variable, suggestion));
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

        % if untreated not selected 
        if ~any(strcmp(chosen_drugs,"Untreated"))
            % I forget which I called it sorry 
            remove_after_TVN = true;
            remove_after_tvn = true;
        else
            remove_after_TVN = false;
            remove_after_tvn = false;
        end 
        
        % want to put on the timecourse data if we're doing that
        if apply_timecourse_data
            load('./workspaces/all_timecourse_doseresponse.mat')
            
            table_combine_temp = struct;
            table_combine_temp.final_data_table = final_data_table;
            
            table_combine_temp.time = timecourse_data;
            
            big_data = column_validation(table_combine_temp);

            final_data_table = vertcat(big_data{:});
            
            choice_variable = unique(final_data_table.(choice_extension))';
        end 

        %prompt the user to select drugs to apply if they want to 
        [indexes_chosen, any_chosen_tf] = listdlg('ListString',choice_variable,'Name',"Select Drugs to apply",'ListSize',[220 380],...
            'CancelString','No Selection'); 

        % Get a list of drugs to apply
        if any_chosen_tf
            drugs_to_apply = choice_variable(indexes_chosen);
            applied_drug = drugs_to_apply{1};
        end

    end 
    %% Verify and initialize drug selections
    filter_drugs

    if apply 
        final_data_table = vertcat(final_data_table,drug_to_apply);
    end 

    %% Remove all img column b/c it buggy
    try
        final_data_table = removevars(final_data_table,"ALL_IMG");
        disp("Removed ALL_IMG column")
    catch
        disp("No ALL_IMG column to remove")
    end

    %% Recuperate our warriors
    merged_table_setup
    drugs = unique(final_data_table.DRUG).';
    ids = unique(final_data_table.ID).';

    %% Sort table by id
    final_data_table = sort_by_column_name(final_data_table, 'ID');

    % need to convert Untreated to remove the 3x and 025x extensions 
    if ~flattened_joint
        % use ID to find untreated inds
        untreated_indicies = find(strcmp(final_data_table.ID, 'Untreated'));
        final_data_table.DRUG(untreated_indicies) = {'Untreated'};

        %%% from here on out, this is basically just individual
        
        % reduce number of untreated to 80 at random
        
        % before getting started, need to randomly select 80 of the
        % untreated
        rng('shuffle')
        unt_rows = find(strcmp(final_data_table.DRUG,"Untreated"));
        unt_data = final_data_table(unt_rows,:);

        % remove rows to add back later 
        final_data_table(unt_rows,:) = [];

        % sample bby
        reduced_unt_rows = datasample(unt_data,num_rows,'Replace',false);

        % now put that thing back where it came from or so help me
        final_data_table = vertcat(final_data_table,reduced_unt_rows);

        disp("Selected 80 sample untreated") 

        
    % if flattened, do dose response analysis to create flattened table 
    else 
        
        %% alright here we need to play around with the untreated for the bayes stuff
        % at this point, Untreated_3x and Untreated_025x are separate entries
        % want to randomly select 80 of each

        rng('shuffle')

        untreated_3x_rows = find(strcmp(final_data_table.DRUG,"Untreated_3x"));
        untreated_025x_rows = find(strcmp(final_data_table.DRUG,"Untreated_025x"));

        untreated_3x_data = final_data_table(untreated_3x_rows,:);
        untreated_025x_data = final_data_table(untreated_025x_rows,:);

        % remove rows from original table to add back later
        final_data_table([untreated_3x_rows; untreated_025x_rows],:) = [];

        % randomly sample some rows
        reduced_3x_rows = datasample(untreated_3x_data,num_rows,'Replace',false);
        reduced_025x_rows = datasample(untreated_025x_data,num_rows,'Replace',false);

        % and put them back in the table
        final_data_table = vertcat(final_data_table,reduced_3x_rows,reduced_025x_rows);

        %% Now flatten the table 
        final_data_table = sort_by_column_name(final_data_table, 'ID');

        % prepare table for flattened table 
        chosen_doses = {'025x', '3x'};
        chosen_ids = ids;

        [doses, doses_by_id, drugs_by_id] = get_dose_list(drugs, chosen_ids); %I'd rewrite this function by cross-referencing rows of get_drug_info()—the way it works now isnt very robust


        % put everything in the flattened table 
        doseResponse_analysis
    end 
    % Run PCA Analysis
    show_pca = false;
    PCA_analysis

else 
    % not flattened joint 
    flattened_joint= false;
    % settings for knn 
    make_median = true;
    make_full = false;
    after_pca = true;
    after_TVN = false;
    % if 3x
    if x3_dose
        chosen_workspace = high_dose_wksp;
        chosen_workspaces = {high_dose_wksp};
        load(strcat('./workspaces/',high_dose_wksp,'.mat'))
        
        % if on later runs and applying condition data, need to make sure
        % to load those workspaces too 
        if runct > 1 && apply_condition_data
            load('./workspaces/all_condition_data_with_unt.mat')
            
            table_combine_temp = struct;
            table_combine_temp.final_data_table = final_data_table;
            if apply_but
                table_combine_temp.but_table = but_table;
            end 
            if apply_chol
                table_combine_temp.chol_table = chol_table;
            end
            if apply_ph
                table_combine_temp.ph_table = ph_table;
            end
            big_data = column_validation(table_combine_temp);

            final_data_table = vertcat(big_data{:});
                    
        end 
    % if not 3x, then 025x
    else
        chosen_workspace = low_dose_wksp;
        chosen_workspaces = {low_dose_wksp};
        load(strcat('./workspaces/',low_dose_wksp,'.mat'))
    end 

    % safety if this doesn't exist
    if ~exist('remove_after_TVN','var')
        remove_after_tvn = true;
        remove_after_TVN = true; %lol i hat emyself 
    end 
    
    % safety declaration 
    choice_extension = "DRUG";
    choice_variable = unique(final_data_table.(choice_extension))';

    %% select our drugs  
    % if we have not already looped and selected our chosen drugs 
    if ~exist('chosen_drugs','var')
        % get a list of all of the drugs 
        all_drugs = unique(final_data_table.(choice_extension))';
        % by default, show what was suggested in bayes_backwards_select
        default_selection = find(ismember(all_drugs, suggestion));
        
        % Prompt user to select a list of drugs
        indexes_chosen = listdlg('ListString',all_drugs,...
            'Name',"Select Drugs",'ListSize',[220 380],'PromptString',"Select Drugs",...
            'InitialValue', default_selection); 

        % Get a list of the drugs chosen
        chosen_drugs = all_drugs(indexes_chosen);

        % if untreated was not selected 
        if ~any(strcmp(chosen_drugs,"Untreated"))
            % I forget which I called it sorry this is the worst
            remove_after_TVN = true;
            remove_after_tvn = true;
        else
            remove_after_TVN = false;
            remove_after_tvn = false;
        end 

        %%% want to add a way to select a different space to apply from
        %%% add in the condition drugs 
        
        if x3_dose && apply_condition_data
            % throw an error if you're missing the workspace 
            if ~isfile('./workspaces/all_condition_data_with_unt.mat')
               error("go download the all_condition_data_with_unt.mat workspace from the box") 
            end
            load('./workspaces/all_condition_data_with_unt.mat')
            variables_to_create = {'apply_but','apply_chol','apply_ph'};
            default_values =      [false,       false,      false];
            target_values = checkboxList("Select Settings", variables_to_create, default_values);

            initialize_variables
            % add this to final data table for each condition, recalculate
            % choice_variable after
            
            table_combine_temp = struct;
            table_combine_temp.final_data_table = final_data_table;
            
            if apply_but
                 table_combine_temp.but_table = but_table;
            end 
            if apply_chol
               table_combine_temp.chol_table = chol_table;
                
            end
            if apply_ph
               table_combine_temp.ph_table = ph_table;
            end 
            
            
            big_data = column_validation(table_combine_temp);

            final_data_table = vertcat(big_data{:});
            
            choice_variable = unique(final_data_table.(choice_extension))';
            
        end 
        
        
        %prompt the user to select drugs to apply if they want to 
        [indexes_chosen, any_chosen_tf] = listdlg('ListString',choice_variable,'Name',"Select Drugs to apply",'ListSize',[220 380],...
            'CancelString','No Selection'); 

        % Get a list of drugs to apply
        if any_chosen_tf
            drugs_to_apply = choice_variable(indexes_chosen);
            applied_drug = drugs_to_apply{1};
            disp(strcat("Applying ", applied_drug))
        end
        
        
    end 

    %% randomly select 80 untreated 
    % before getting started, need to randomly select 80 of the
    % untreated
    rng('shuffle')
    unt_rows = find(strcmp(final_data_table.DRUG,"Untreated"));
    unt_data = final_data_table(unt_rows,:);

    % remove rows to add back later 
    final_data_table(unt_rows,:) = [];

    % sample bby
    reduced_unt_rows = datasample(unt_data,num_rows,'Replace',false);

    % now put that thing back where it came from or so help me
    final_data_table = vertcat(final_data_table,reduced_unt_rows);

    disp("Selected 80 sample untreated") 
    
    % do PCA analysis to filter drugs and set up workspace for backwards
    % selecting 
    PCA_analysis



end 
%% make the bayes struct only if this is the first time around 
% if it's the first run...
if runct == 1
    % if joint profile with flattened table 
    if do_joint_profile && flattened_joint
        
        % don't put untreated in the struct if we're gonna remove it
        % later
        if remove_after_TVN
            disp("after initial analysis, remove_after_TVN true")            
            all_vals = unique(final_data_table.ID);
            disp("value of all_vals")
            disp(all_vals)
            % find the untreated
            unt_loc = find(strcmp(all_vals,"Untreated"));
            all_vals(unt_loc) = [];
            bayes_struct = make_bayes_struct(all_vals);
        else 
            bayes_struct = make_bayes_struct(unique(final_data_table.ID));
        end 
        
        
    % if individual profile 
    else
        % don't put untreated in the struct if we're gonna remove it
        % later
        if remove_after_TVN
            disp("after initial analysis, remove_after_TVN still true")
            all_vals = unique(final_data_table.DRUG);
            disp("value of all_vals")
            disp(all_vals)
            % find the untreated
            unt_loc = find(strcmp(all_vals,"Untreated"));
            all_vals(unt_loc) = [];
            bayes_struct = make_bayes_struct(all_vals);
        else
            bayes_struct = make_bayes_struct(unique(final_data_table.DRUG));
        end 
        
    end
    
    % make something to keep track of the variables
    all_variable_sets = cell(maxrun,1); 
    %% if doing random replicates, this is where we swap labels 
    % regardless if individual or joint, this will work
    % on this first run is where we want to determine our random labels
    if replicates_random 
        rng('default')
        % don't want to do untreated
        % can skip over it both in the flat and not flat version this way 
         unt_rows = strcmp(final_data_table.DRUG,"Untreated") | strcmp(final_data_table.ID,"Untreated");
         untreated_table = final_data_table(unt_rows,:);
         no_untreated_final_data = final_data_table(~unt_rows,:);
         % get random labels
         % if we are applying a drug, don't want to swap the labels for
         % that 
         if apply
             % same process, find inds of applied drug and remove it so it
             % doesn't get swapped and then put that thing back where it
             % came from 
             applied_drug_inds = strcmp(no_untreated_final_data.DRUG,applied_drug) | strcmp(no_untreated_final_data.ID,applied_drug);
             applied_table = no_untreated_final_data(applied_drug_inds,:);
             no_unt_or_applied = no_untreated_final_data(~applied_drug_inds,:);
             
             % this is the random order that we need to use and apply each time!!!
             REPLICATE_RANDOM_ORDER = randperm(length(no_unt_or_applied.DRUG));
             % reassign 
             no_unt_or_applied.DRUG = no_unt_or_applied.DRUG(REPLICATE_RANDOM_ORDER,:);
             
             % now put the final data table back together

             final_data_table = vertcat(no_unt_or_applied, untreated_table,applied_table);
             
         else  
         
             % this is the random order that we need to use and apply each time!!!
             REPLICATE_RANDOM_ORDER = randperm(length(no_untreated_final_data.DRUG));
             % reassign 
             no_untreated_final_data.DRUG = no_untreated_final_data.DRUG(REPLICATE_RANDOM_ORDER,:);

             % now put the final data table back together

             final_data_table = vertcat(no_untreated_final_data, untreated_table);
         end 
    end 
else
%% if it's not the first run and we're doing random, apply the random labels here
    if replicates_random 
        unt_rows = strcmp(final_data_table.DRUG,"Untreated") | strcmp(final_data_table.ID,"Untreated");
        untreated_table = final_data_table(unt_rows,:);
        no_untreated_final_data = final_data_table(~unt_rows,:);
        
        if apply
            % again making sure that the applie drug does not get its
            % labels changed 
            applied_drug_inds = strcmp(no_untreated_final_data.DRUG,applied_drug) | strcmp(no_untreated_final_data.ID,applied_drug);
            applied_table = no_untreated_final_data(applied_drug_inds,:);
            no_unt_or_applied = no_untreated_final_data(~applied_drug_inds,:);
            no_unt_or_applied.DRUG = no_unt_or_applied.DRUG(REPLICATE_RANDOM_ORDER,:);
              
            % now put the final data table back together
            final_data_table = vertcat(no_unt_or_applied, untreated_table,applied_table);
            
        else
            % apply the random order 
            no_untreated_final_data.DRUG = no_untreated_final_data.DRUG(REPLICATE_RANDOM_ORDER,:);
            final_data_table = vertcat(no_untreated_final_data, untreated_table); 
        end    
    end 
end 

% set knn settings -- don't create any of the graphs
create_cm = 'off';
show_pca = false;
create_graph = false;

%% now do knn analysis to finish setting up workspace 

KNN_analysis
%% run backwards select w/o overhead gui and auto save at end 

disp("beginning backwards_select_analysis")

% run through backwards select with the currently set up workspace 
backwards_select_analysis

%% Pull out best result from backwards select 
% after backwards_select_analysis, need to pull out the best result
% from the overall vars --> just straight up look for highest percent
% and if there is more than one, take lowest value

% get all of the fieldnames for overall vars 
field_names = fieldnames(overall_vars);

% only find variable set fields (ignore DRUGS and other metadata)
is_var = contains(field_names,'vars');

% only get var fields
var_fields = field_names(is_var);

% all are in the form vars_18_pct_40
% jut get the pct numerical values
all_ov_pcts = str2double((extractAfter(var_fields,"pct_")));

% now want to find the max value in all_ov_pcts
[max_ov_value] = max(all_ov_pcts);

% want to make sure we get all of the possible max indicies
max_ov_idx = find(all_ov_pcts == max_ov_value);

% grab the one with the fewest variables, so the one closest to the
% bottom
if length(max_ov_idx) > 1
    max_ov_idx = max_ov_idx(end); 
end 

% get the best ov field 
ov_field = var_fields{max_ov_idx};
disp(ov_field)
% get corresponding feature set
starting_vars = overall_vars.(ov_field);

% save in all_variable_sets
all_variable_sets{runct} = starting_vars;

% from here, need to do analysis again with this set and save the knn
% results 
do_feature_reduction = false;

%% run the analysis again back with the best set

% perform TVN with newly selected features
% but also don't do that if we're not doing that 
if no_TVN
    pca_whitened_table = final_data_table; 
else
    TVN_transform
end 

%  appropriately remove untreated if needbe 
if remove_after_TVN
    disp("after selection of optimal variables, remove_after_TVN is true")
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
    disp("removed unt")
    disp(drugs)
end 

%% if applied drug, need to remove here 

if apply
    whitened_table_copy = pca_whitened_table;         
    for qi = addedDrug
        added_drug_1 = qi{1};
        applied_inds = find(strcmp(pca_whitened_table.(choice_extension),added_drug_1));
    % remove!!!
        pca_whitened_table(applied_inds,:) = [];
    end 

end 

%% start PCA calcs

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,pca_whitened_table,'OutputFormat', 'uniform');

%find column names
numeric_final_col_names = pca_whitened_table.Properties.VariableNames(numeric_final_data_cols);

% get numeric data
ndata = pca_whitened_table{:,numeric_final_data_cols};

% non numeric data
non_numeric = pca_whitened_table(:,~numeric_final_data_cols);


%% perform PCA calculation 
[coeff,scores,pcvars] = pca(ndata);
vari = cumsum(pcvars)./sum(pcvars);
coeff=coeff';
pcvari = pcvars/sum(pcvars);

% put data in a table for better organization
scores_table = horzcat(non_numeric,array2table(scores)); 

% if we're applying, all we care about is the applied version
if apply
    % use whitened_table_copy

    ndata2 = whitened_table_copy{:,numeric_final_data_cols};  

    non_numeric = whitened_table_copy(:,~numeric_final_data_cols);
    
    %apply eigevnectors
    newspace = ndata2*coeff';
    
    % organize newspace into a table
    scores_table = horzcat(non_numeric,array2table(newspace));

    
end 

% get list of drugs 
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
    if ~flattened_joint && do_joint_profile
        [drug_medians,drugs] = table_median(scores_table,'DRUG');
    else 
        [drug_medians,drugs] = table_median(scores_table,choice_extension);
    end 
else
    % if not, want to use pca_whitened_table
    [drug_medians,drugs] = table_median(pca_whitened_table,choice_extension);
end

%%% random here

if random_at_end
    % first time need to set up the random labels 
    if runct == 1
        END_RANDOM_ORDER = randperm(length(drugs));
        drugs = drugs(END_RANDOM_ORDER);
    else
        % just apply the random permutation already determined 
        drugs = drugs(END_RANDOM_ORDER);
        
    end 
    
end 

% use medians as knn_data
knn_data = drug_medians;
knn_drugs = drugs';

% save_bayes_results
save_bayes_results = true;

% give some user update
disp("running with")
disp(ov_field)

% perform knn, saving into bayes struct 
knn_helper

disp(pct_correct)