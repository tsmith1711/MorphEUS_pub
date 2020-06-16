%% Kruskal-Wallis test Multicompare Analysis with user prompts for workspace, and label selection

% Performs KW significance test between selected data sets. Can be done
% with both population data and individual cell data
% can also create boxplots showing data distribution, excluding outliers
% futher options are explained in Prompt Settings 


%% Add Paths
add_all_paths

% turn off latex interpeter
set(groot, 'defaultAxesTickLabelInterpreter', 'none')

figure_folder_path =  "./figures/";

% antiquated: removed a set of poorly imaged replicates from a previous
% data set
remove_bad_treated = false;
% antiquated: removed extra controls for ease of analysis, no longer part o
% the process
remove_extra_controls = false;

%% Prompt settings

variables_to_create = {'by_cell','large_group', 'remove_INH_control', 'efflux_in_lipid', 'violin_plot','rearrange_mode','make_boxplot','conditions_defaults'};
default_values =      [true,     true,          true,                  true,               false         false,         true,            false];
disp("VARIABLE EXPLANATION")
disp("--")
disp("by_cell = perform statistical analysis on workspaces in cell_workspaces.")
disp("large_group = use the broader categories for drug pathways.")
disp("remove_INH_control = automatically remove INH_control from the list.")
disp("efflux_in_lipid = move drugs previously classified as efflux into the lipid category.")
disp("violin_plot = create violin plots")
disp("rearrange_mode = allow user to rearrange the order of the drugs in the plot")
disp("make_boxplot = create boxplots")
disp("conditions_defaults = pre-select defaults for creating the condition comparison plots")
disp("--")

target_values = checkboxList("Select Multicompare Settings", variables_to_create, default_values);

initialize_variables

%% make this global to make it accessible inside ver_figsave()
global chosen_workspaces


%% Find list of workspaces, prompt user for choice, and load chosen

% if user selected by cell, want to look in cell_workspaces directory
if by_cell
    workspace_directory = './cell_workspaces';
else
    workspace_directory = './workspaces';
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
    'Name',"Select workspace",'ListSize',[340 380],...
    'InitialValue', default_selection);

% if the user did not select anything, exit (there is a catch for this in load_workspaces too)
if ~any_chosen_tf
    return
end

% use the index to get the list of user-selected workspaces
chosen_workspaces = workspace_list(indexes_chosen);

%% run the script to load the workspaces
load_workspaces

merged_table_setup

% assign moas based on user selection of large group earlier (will only
% work for 52 drug set right now, could maybe make it so that it adds a new column for MoA?)

% if not cell mode, assign moas (would take too long in cell mode)
if ~by_cell
    final_data_table = assign_moas(final_data_table,large_group);
end 

final_data_table = sortrows(final_data_table,2);

%% ask if user would like to compare by drug or MoA or whatevr

ways_to_plot = non_numeric.Properties.VariableNames;

% prompt the user to select what they would like to group by
[indexes_chosen, any_chosen_tf] = listdlg('ListString',ways_to_plot,...
    'Name',"Select grouping method",'ListSize',[220 380],'SelectionMode','single');

% if the user did not select anything, exit (there is a catch for this in load_workspaces too)
if ~any_chosen_tf
    return
end

chosen_method = ways_to_plot(indexes_chosen);


% get the table variable to use 
chosen_method = chosen_method{1};

% find the indexes for the chosen method
chosen_values = unique(final_data_table.(chosen_method))';

chosen_method_indexes = cell2struct(cell(1,length(chosen_values)), chosen_values, 2);

for i = chosen_values
    drug = i{1};
    chosen_method_indexes.(drug) = find(strcmp(final_data_table.(chosen_method), drug)).';
end


%% From the list of drugs, prompt user to select which drugs to use

if conditions_defaults
    correct_order =  {'Untreated_x3','Untreated_but','Untreated_chol', 'Untreated_ph',...
    'EMB_3x_x3','EMB_but','EMB_chol','EMB_ph','INH_3x_x3','INH_but','INH_chol',...
    'INH_ph', 'Pre_3x_x3','Pre_but','Pre_chol','Pre_ph','BDQ_3x_x3','BDQ_but',...
    'BDQ_chol','BDQ_ph','Clz_3x_x3','Clz_but','Clz_chol','Clz_ph','Mox_3x_x3',...
    'Mox_but','Mox_chol','Mox_ph','Lin_3x_x3','Lin_but','Lin_chol','Lin_ph',...
    'RIF_3x_x3','RIF_but','RIF_chol','RIF_ph','RifT_3x_x3','RifT_but','RifT_chol',...
    'RifT_ph','PZA_3x_x3','PZA_but','PZA_chol','PZA_ph'};
    default_selection = find(ismember(chosen_values, correct_order));
    
else 
    
    default_selection = [];
end 
% Prompt user to select a list of drugs
[indexes_chosen, any_chosen_tf] = listdlg('ListString',chosen_values,...
    'Name',"Select which to compare",'ListSize',[220 380],'PromptString',"Select Drugs",...
    'InitialValue',default_selection); 

% if the user did not select anything, exit
if ~any_chosen_tf
    return
end

% Get a list of the drugs chosen
chosen_to_analyze = chosen_values(indexes_chosen);

%% remove non-selected values

% drugs user decided to ignore
removed = chosen_values;
indx = ismember(chosen_values, chosen_to_analyze);
removed(indx) = [];
clear indx;

% remove ignored drugs from the list 
removallist = [];
for i = removed
    drug = i{1};
    removallist = [chosen_method_indexes.(drug) removallist];
end

final_data_table(removallist,:) = [];

% find the indexes for the chosen method
chosen_values = unique(final_data_table.(chosen_method))';

chosen_method_indexes = cell2struct(cell(1,length(chosen_values)), chosen_values, 2);

for i = chosen_values
    drug = i{1};
    chosen_method_indexes.(drug) = find(strcmp(final_data_table.(chosen_method), drug)).';
end

% include way to get cell count from this table 
if by_cell
    % use a syto feature and an fm feature and a shape feature
    % we could really just input static values but idk I want it to be
    % dynamic b/c why not 
    bact_feature = contains(final_data_table.Properties.VariableNames,"FEATURE");
    % keep the true columns
    bact_feature = find(bact_feature >0);
    % we only need one, so just grab the first 
    bact_feature = bact_feature(1);
    
    
    fm_feature = contains(final_data_table.Properties.VariableNames,"f_");
    % keep the true columns
    fm_feature = find(fm_feature >0);
    % we only need one, so just grab the first 
    fm_feature = fm_feature(1);
    
    syto_feature = contains(final_data_table.Properties.VariableNames,"s_");
    % keep the true columns
    syto_feature = find(syto_feature >0);
    % we only need one, so just grab the first 
    syto_feature = syto_feature(1);
    
    %% now use those features to get a cell count (only want nonzero values for syto/fm)
    % need to loop through each 
    
    bact_cell_counts = struct;
    FM_cell_counts = struct;
    syto_cell_counts = struct;
   
    for w = chosen_to_analyze
        drug = w{1};
        bact_size = size(final_data_table(chosen_method_indexes.(drug),bact_feature),1);
        
        % need nonzero values
        FM_size = size(find(final_data_table{chosen_method_indexes.(drug),fm_feature} > 0),1);
        
        % need nonzero values
        syto_size = size(find(final_data_table{chosen_method_indexes.(drug),syto_feature} > 0),1);
        
        
        % save all values in structure
        bact_cell_counts.(drug) = bact_size;
        FM_cell_counts.(drug) = FM_size;
        syto_cell_counts.(drug) = syto_size;
    
    end 
    
    
end 

%% if we're doing condition data, let's go ahead and rearrange the way we want b/c it's hard later
if conditions_defaults


    correct_order =  {'Untreated_x3','Untreated_but','Untreated_chol', 'Untreated_ph',...
    'EMB_3x_x3','EMB_but','EMB_chol','EMB_ph','INH_3x_x3','INH_but','INH_chol',...
    'INH_ph', 'Pre_3x_x3','Pre_but','Pre_chol','Pre_ph','BDQ_3x_x3','BDQ_but',...
    'BDQ_chol','BDQ_ph','Clz_3x_x3','Clz_but','Clz_chol','Clz_ph','Mox_3x_x3',...
    'Mox_but','Mox_chol','Mox_ph','Lin_3x_x3','Lin_but','Lin_chol','Lin_ph',...
    'RIF_3x_x3','RIF_but','RIF_chol','RIF_ph','RifT_3x_x3','RifT_but','RifT_chol',...
    'RifT_ph','PZA_3x_x3','PZA_but','PZA_chol','PZA_ph'};

    % oops
    correct_order = flip(correct_order);

    copy_table = final_data_table;
    % empty version of final_data_table to fill in with rearranged version 
    copy_table(:,:) = [];
    for o = correct_order
        current_dexp = o{1};
        % find all rows with that values
        
        inds = find(strcmp(final_data_table.DRUG_EXP,current_dexp));
        
        copy_table = vertcat(copy_table,final_data_table(inds,:));
        
    end 
    
    final_data_table = copy_table;
end 



%% allow user to select variables
with_percent = [numeric_final_col_names {'percent_FM'} {'percent_syto'}];

% Prompt user to select features 
[indexes_chosen, any_chosen_tf] = listdlg('ListString',with_percent,...
    'Name',"Select feature",'ListSize',[220 380],'PromptString',"Select Feature"); 

% if the user did not select anything, exit
if ~any_chosen_tf
    return
end
clear any_chosen_tf

chosen_comparison_variable = with_percent(indexes_chosen);

first_loop = true;


for j = chosen_comparison_variable
    % loop through each var
    chosen_var = j{1};
    % and multicompare it bby
    % get the labels
    final_label_selection = final_data_table.(chosen_method);
    % get the data 
    if ~strcmp(chosen_var,"percent_FM") && ~strcmp(chosen_var,"percent_syto")
        final_data_selection = final_data_table.(chosen_var);
        
        % need to convert to percent if that's what we want 
    elseif strcmp(chosen_var,"percent_FM")
        final_data_selection = final_data_table.FEATURE_2_count > 0;
        final_data_selection = double(final_data_selection);
        % really we want to just get that percent value 
        % loop through all different drugs and output the percent, maybe?
        all_labels = unique(final_label_selection,'stable')';
        
        FM_percents = zeros(1,length(all_labels));
        ct = 1;
        for k = all_labels
            current_drug_label = k{1};
            inds = find(strcmp(final_label_selection,current_drug_label));
            % get the percent stained 
            pct_val = mean(final_data_selection(inds))*100;
            FM_percents(ct) = pct_val;
            ct = ct+1;
        end 
        
        FM_stain_table = table(all_labels,FM_percents,'VariableNames',{'DRUG','PCT_STAIN'});
    
    elseif strcmp(chosen_var,"percent_syto")
        final_data_selection = final_data_table.FEATURE_1_count > 0;
        final_data_selection = double(final_data_selection);
        % really we want to just get that percent value 
        % loop through all different drugs and output the percent, maybe?
        all_labels = unique(final_label_selection,'stable')';
        
        syto_percents = zeros(1,length(all_labels));
        ct = 1;
        for k = all_labels
            current_drug_label = k{1};
            inds = find(strcmp(final_label_selection,current_drug_label));
            % get the percent stained 
            pct_val = mean(final_data_selection(inds))*100;
            syto_percents(ct) = pct_val;
            ct = ct+1;
        end 
        
        syto_stain_table = table(all_labels,syto_percents,'VariableNames',{'DRUG','PCT_STAIN'});
    
        
    end 
    
    
    % need to do something with the s and f data for the by cell data since
    % there are a lot of extra zeros
    if by_cell
        if contains(chosen_var,"s_") || contains(chosen_var,"f_")
            % find inds of nonzero values
            nonzero = find(final_data_selection > 0);
            
            % only keep nonzero values
            final_data_selection = final_data_selection(nonzero);
            % do same for labels so it still matches
            final_label_selection = final_label_selection(nonzero);
            
        end 
    end
  %  [p,t,stats] = anova1(final_data_selection,final_label_selection,'off');
    [p,t,stats] = kruskalwallis(final_data_selection,final_label_selection,'off');
    
    if make_boxplot
        
        % don't want boxplot for percetn FM 
        if strcmp(chosen_var,"percent_FM")
            figure
            scatter(1:length(FM_stain_table.PCT_STAIN),FM_stain_table.PCT_STAIN,'filled')
            xticks(1:length(FM_stain_table.DRUG))
            xticklabels(FM_stain_table.DRUG)
            view([90 -90])
            % or percent syto
        elseif strcmp(chosen_var,"percent_syto")
            figure
            scatter(1:length(syto_stain_table.PCT_STAIN),syto_stain_table.PCT_STAIN,'filled')
            xticks(1:length(syto_stain_table.DRUG))
            xticklabels(syto_stain_table.DRUG)
            view([90 -90])
        else 
            figure
            boxplot(final_data_selection,final_label_selection,'BoxStyle','filled',...
                'Colors',[0.2 0.2 0.8],'Symbol','')

             lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
            set(lines, 'Color', [0.95 0.1 0.1]);
            pause(0.01)
            view([90 -90])
        end 
        title(chosen_var,'Interpreter','none')
    end 

    disp(strcat("RESULTS FROM ONE WAY ANOVA FOR ", chosen_var))
    
    % want to give user chance to rearrange their multicompare if they want
    % to
    if rearrange_mode && first_loop
        
   
        % make a list of numbers for each name
        numbers = 1:length(stats.gnames);

        % correct the orientation 
        numbers = numbers';

        % concatenate the numbers with the names and display
        number_names = cellstr(strcat(num2str(numbers), " ", stats.gnames));
        
        % need to display it in figure so we can see all options
        copy_stats = stats;
        
        % use the number indicies 
        copy_stats.gnames = number_names;
        figure
        multcompare(copy_stats);
        
        % prompt user to write new order 
        answer = inputdlg('Enter space-separated new order:',...
                     'Sample', [1 50]);
                 
        % if user cancelled, set rearrange mode to false and move on
        % but if they did not cancel, rearrange!
        if ~isempty(answer)
            user_val = str2num(answer{1});

            rearrangeMulticompare
            title({chosen_var, chosen_workspace},'Interpreter','none')
        else 
            rearrange_mode = false;
            figure 
            [c,m,h,nms] = multcompare(stats);
            title({chosen_var, chosen_workspace},'Interpreter','none')
        end 
        
        % if we are on a later loop 
    elseif rearrange_mode && ~first_loop
        rearrangeMulticompare
        title({chosen_var, chosen_workspace},'Interpreter','none')

        % if we're not in rearrange mode, just do the normal thing 
    else
        % create a new figure each time
        figure 
        [c,m,h,nms] = multcompare(stats);
        title({chosen_var, chosen_workspace},'Interpreter','none')
        
%         figure 
%         [c2,m2,h2,nms2] = multcompare(stats2);
%         title({chosen_var, chosen_workspace},'Interpreter','none')
    end 
    

    %ver_figsave(h,strcat(figure_folder_path,chosen_comparison_variable, " multicompare "))

    % want to convert the c array to a more readable table format

    comparison_table = multicompare_table(c,nms)
    
    
    first_loop = false; 
end 

% if cell mode, want to display the cell counts we calculated
if by_cell
    bact_cell_counts
    FM_cell_counts
    syto_cell_counts 
    
end

%% if we're doing violin plots 
if by_cell && violin_plot
    for j = chosen_comparison_variable
        chosen_var = j{1};
        data_array = {};
        % need to loop through each drug: data has to be set up differently
        % in this one
        all_categories = unique(final_data_table.(chosen_method))';
         
        final_data_selection = final_data_table.(chosen_var);
        final_label_selection = final_data_table.(chosen_method);
        if by_cell
            if contains(chosen_var,"s_") || contains(chosen_var,"f_")
                % find inds of nonzero values
                nonzero = find(final_data_selection > 0);

                % only keep nonzero values
                final_data_selection = final_data_selection(nonzero);
                % do same for labels so it still matches
                final_label_selection = final_label_selection(nonzero);

            end 
        end
        
        ct = 1;
        for w = all_categories
            category = w{1};
            
            % get the rows for this category
            correct_inds = find(strcmp(final_label_selection,category));
            
            % get those rows for that selected variable
            the_data = final_data_selection(correct_inds);
            
            data_array(ct) = {the_data};
            ct = ct+1;
        end 
        figure
        [h,L,MX,MED,bw]=violin(data_array);
        xticklabels(strrep(all_categories,"_"," "))
        xtickangle(45)
        title(chosen_var,'Interpreter','none')
    end 
    
end 

