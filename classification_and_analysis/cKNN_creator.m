%% Create a cKNN for a bayes workspace
%%% REQUIRES: existing bayes structure created with classification_trials
% by default, searches trial_results, allows user to make selection
% will create heatmaps and connectivity maps
% if selecting crossval_mode, will not create heatmaps/connectivity maps,
% will just read out what the highest MoA and nearest neighbor for each
% drug is 

%clear
% turn off tex interpreter for axes 
set(groot, 'defaultAxesTickLabelInterpreter', 'none')

%% let user load workspace of choice

% add obligatory paths
add_all_paths

disp("Select a bayes workspace")
[bayes_file,path] = uigetfile('./trial_results/*.mat');
% if user did not make a selection, exit
if isequal(bayes_file,0)
   disp('User selected Cancel');
   return
end

% load selected workspace
load(fullfile(path,bayes_file))
disp(['Loading ', fullfile(path,bayes_file)]);


%% new map - bw 
vec = [100; 0;];
raw = [ .9 .9 .9; 0 0 0];
N = 128;
bw_map = interp1(vec,raw,linspace(100,0,N),'pchip');

%% new map - red 
vec = [100; 50; 0;];
raw = [ 1 1 1; ...
    255/255 104/255 104/255; ...
    47/255 0 0];
N = 128;
red_map = interp1(vec,raw,linspace(100,0,N),'pchip');


%% new map -> parula plus some purple 

parula_plus = parula;

parula_plus = parula_plus(6:end,:);

% need to add purple color to beginning -> fade from first blue color to a
% purple color
first_blue = parula_plus(1,:);
purple_hue = [59/255 14/255 101/255];
vec = [0;100];
raw = [first_blue; purple_hue];
N = 14; 
blue_to_purple = interp1(vec,raw,linspace(100,0,N),'pchip');

parula_plus = [blue_to_purple; parula_plus];
% get a purple hue to fade to
parula2 = parula_plus(10:end,:);

%% initialize variables 
variables_to_create = {'large_group','use_drug_full_name','crossval_mode'};
default_values =      [true,          true,              false];
disp("VARIABLE EXPLANATION")
disp("--")
disp("large_group = use the broader categories for drug pathways.")
disp("use_drug_full_name = use full names instead of 3-4 letter abbreviations")
disp("crossval_mode = don't display figures, output a list of highest connections for bidirectional")
disp("--")
target_values = checkboxList("Select PCA_analysis Settings", variables_to_create, default_values);

initialize_variables %Initialize variables_to_create as target_values

% sets up color scheme 
red_colormap = true;
chosen_colormap = red_map;


efflux_in_lipid = true;

% external script to
% create the colormap, all_categories, and color_struct
knn_color_setup

%% convert bayes structure to a percent table

% creates a one-way table
table1 = struct2table(bayes_struct);

data_as_array = [];

% go throughe each column
for i = 1:size(table1,2)
    % current
    current  = table1{:,i};
    new_col = struct2array(current)';
    total_runs = sum(new_col);
    % re scale each as % of total 
    new_col_pct = new_col/total_runs*100;
    
    data_as_array = horzcat(data_as_array,new_col_pct);
    
end 

run_tot_string = strcat(num2str(total_runs), " total runs");
disp(run_tot_string)


%% convert into table form with row and column names 

bayes_table_pct = array2table(data_as_array,'VariableNames',table1.Properties.VariableNames,'RowNames',table1.Properties.VariableNames);

% change to all caps or full names if we're into that 

if use_drug_full_name
    bayes_table_pct.Properties.VariableNames = replace_full_name(bayes_table_pct.Properties.VariableNames);
else
    
    bayes_table_pct.Properties.VariableNames = replace_allcap_name(bayes_table_pct.Properties.VariableNames);
end 

%% need to restructure for digraph: loop through each

all_drugs = bayes_table_pct.Properties.VariableNames;

starting_drugs = [];
neighbor_drugs = [];
weights = [];

% loop through all drugs
for i = all_drugs
    drug = i{1};
    
    drug_col = bayes_table_pct.(drug);
    
    % want to find nonzero values
    nonzero_ind = find(drug_col >0);
    nonzero_vals = drug_col(nonzero_ind);
    neighbors = all_drugs(nonzero_ind)';
    
    % need to make a drug col of the corresponing length where all values
    % are just the starting drug
    drugs = neighbors;  
    for j = 1:length(drugs)
        drugs{j} = drug;
    end 
    
    % add all starting drugs to full array 
    starting_drugs = [starting_drugs; drugs];
    % add all neighbor drugs to full array
    neighbor_drugs = [neighbor_drugs; neighbors];
    % add all weights to full array
    weights = [weights; nonzero_vals];
    
end 


%% put together into a results table 
results_table = table(starting_drugs,neighbor_drugs,weights,'VariableNames', ...
    {'DRUG','NEIGHBOR','DISTANCE'});

% do case insensitive sort 
[~,idx]=sort(upper(results_table.DRUG));
results_table=results_table(idx,:);

%% some colormap settings 
% whether to have a white background or a grey background 
if red_colormap
    wbg = false;
else
    wbg = true;
end 

%% Create interactive connectivity plot with a *~*slider*~*

%connectivity_plot(results_table,color_struct,all_categories,chosen_colormap,wbg)
%sgtitle({bayes_file,run_tot_string},'Interpreter','none')

%% find all of the MoA buddies

results_table.NEIGHBOR_MOA = change_to_moa(results_table.NEIGHBOR,all_categories);

drugs = unique(results_table.DRUG)';

% create empty moa table 
moa_table  = cell2table(cell(0,3), 'VariableNames', {'DRUG', 'NEIGHBOR', 'DISTANCE'});

% loop through each drug 
for i = drugs
    drug = i{1};
    % want to skip untreated and water
    
    % get just rows for that drug
    row_inds = find(strcmp(results_table.DRUG,drug));
    
    drug_rows = results_table(row_inds,:);
    
    % now want to go through each unique MoA, weighting with the weights
    % that are weight
    moas = unique(drug_rows.NEIGHBOR_MOA)';
    
    for j = moas
        moa = j{1};
        
        moa_row_inds = find(strcmp(drug_rows.NEIGHBOR_MOA,moa));
        moa_rows = drug_rows(moa_row_inds,:);
        
        % get sum of weights
        sum_of_weights = sum(moa_rows.DISTANCE);
        
        % make a new row for the table bby
        new_row = table({drug},{moa},sum_of_weights, 'VariableNames', {'DRUG', 'NEIGHBOR', 'DISTANCE'});
        
        moa_table =vertcat(moa_table,new_row);
    end 
    
end 

moa_table = sortrows(moa_table,'DISTANCE','descend');

moa_table = sortrows(moa_table,'DRUG');

%disp(moa_table)

%% make moa table and connectivity plot 

%connectivity_plot(moa_table,color_struct,all_categories,chosen_colormap,wbg)
%sgtitle({bayes_file,run_tot_string},'Interpreter','none')

%% want to go through and make it into funky two way table time

% add a moa column
results_table.MOA = results_table.DRUG;
results_table.MOA = change_to_moa(results_table.DRUG,all_categories);

% resort by moa
results_table = sortrows(results_table,5);

% add moa column for sorting purposes 
moa_table.MOA = moa_table.DRUG;
moa_table.MOA = change_to_moa(moa_table.DRUG,all_categories);

% do case insensitive sort 
[~,idx]=sort(upper(moa_table.DRUG));
moa_table=moa_table(idx,:);
% group by MOA
moa_table = sortrows(moa_table,'MOA');

% use make two way table to create two way table from moa table
new_two = make_two_way_table(moa_table,'DRUG','NEIGHBOR');

%make_connectivity_heatmap(new_two,true);
%colormap(parula_plus)

%% do this by drug as well
% generate two way table 
two_way_table= make_two_way_table(results_table,'DRUG','NEIGHBOR');
% make heatmap 
%make_connectivity_heatmap(two_way_table,false);
%colormap(parula_plus)

%xtickangle(90)

%title("drugs on bottom have the neighbors of drugs on left")

%% want to create a table that adds together all of the connections 
% would it be easier to start with the bayes struct or 1 way table? 
% we'd need to add together first again before taking the percents... hm

% could actually just get raw values back by multiplying and dividng by 100

raw_two_way_table = two_way_table; 

raw_two_way_table{:,:} = raw_two_way_table{:,:}/100*total_runs;

converted_to_half = raw_two_way_table;

converted_to_half{:,:} = converted_to_half{:,:}/(total_runs*2)*100;

% go through each column and add its value to the row version

all_drugs = converted_to_half.Properties.VariableNames;

for i = all_drugs
    drug = i{1};
   
    % now go through all of the rows in just that column
    for j = all_drugs
        current_row_drug = j{1};
        
        current_value = converted_to_half{current_row_drug,drug};
        
        % want to add it to its transpose
        converted_to_half{drug,current_row_drug} = current_value + converted_to_half{drug,current_row_drug};
        
        % and then set the value to zero
        converted_to_half{current_row_drug,drug} = 0;
        
    end 
    
end 

%% make it into a heatmap 
% create heatmap
if ~crossval_mode
    make_connectivity_heatmap(converted_to_half,false);
    xtickangle(90)
    title({bayes_file,"as percent - bidirectional"},'Interpreter','none')
    colormap(parula_plus)
end 
%% create it as the full (not half version)
% this will allow us to do the MoA thing from before
% need to take into account 

full_raw_table = converted_to_half;

symmetric_table = converted_to_half;

full_raw_table{:,:} = converted_to_half{:,:}*(total_runs*2)/100;

bottom_half = tril(full_raw_table{:,:});

bottom_half_pct = tril(converted_to_half{:,:});

symmetric_table{:,:} = symmetric_table{:,:}';

full_raw_table{:,:} = full_raw_table{:,:}';

% now we have this as a full ttable 
full_raw_table{:,:} = full_raw_table{:,:} + bottom_half;

symmetric_table{:,:} = symmetric_table{:,:} + bottom_half_pct;

% row and column totals will be the same
row_totals = sum(full_raw_table{:,:});

all_drugs = full_raw_table.Properties.VariableNames;

% go through all drugs, keeping in mind the number of connections in each
% category
bidirectional_table  = cell2table(cell(0,3), 'VariableNames', {'DRUG', 'NEIGHBOR', 'CONNECTIONS'});
% let's convert to a one way table first 
for i = 1:length(all_drugs)
    drug = all_drugs{i};
    all_connections = row_totals(i);
    
    for j = 1:length(all_drugs)
        drug2 = all_drugs{j};
        % now go through each of its matches
        num_connect = full_raw_table{drug,drug2};
        
        new_row = table({drug},{drug2},num_connect, 'VariableNames', {'DRUG', 'NEIGHBOR', 'CONNECTIONS'});
        
        bidirectional_table = vertcat(bidirectional_table,new_row);
    end 
    
end 


%% make it into a heatmap 
% create heatmap
if ~crossval_mode
    make_connectivity_heatmap(symmetric_table,false);
    xtickangle(90)
    title({bayes_file,"symmetric as percent - bidirectional"},'Interpreter','none')
    colormap(parula_plus)
end 

%% try this


converted_to_half_raw = converted_to_half;

converted_to_half_raw{:,:} = converted_to_half_raw{:,:}*(total_runs*2)/100;

% go through all drugs, keeping in mind the number of connections in each
% category
bidirectional_table2  = cell2table(cell(0,3), 'VariableNames', {'DRUG', 'NEIGHBOR', 'DISTANCE'});
% let's convert to a one way table first 
for i = 1:length(all_drugs)
    drug = all_drugs{i};
    
    for j = 1:length(all_drugs)
        drug2 = all_drugs{j};
        % now go through each of its matches
        num_connect = converted_to_half_raw{drug,drug2};
        
        new_row = table({drug},{drug2},num_connect, 'VariableNames', {'DRUG', 'NEIGHBOR', 'DISTANCE'});
        
        bidirectional_table2 = vertcat(bidirectional_table2,new_row);
    end 
    
end 

% remove rows with zero

zero_inds = find(bidirectional_table2.DISTANCE ==0);

bidirectional_table2(zero_inds,:) = [];

bidirectional_table2.DISTANCE = bidirectional_table2.DISTANCE*100/(total_runs*2);

if ~crossval_mode
    connectivity_plot(bidirectional_table2,color_struct,all_categories,chosen_colormap,wbg)
    sgtitle({bayes_file,strcat(run_tot_string,"-bidirectional")},'Interpreter','none')
end 

%% pause -- clustergram time

% pct_raw_table = full_raw_table;
% pct_raw_table{:,:} = 100*pct_raw_table{:,:}/140; 
% 
% ztest = clustergram(converted_to_half{:,:},'Colormap',parula_plus,...
%     'DisplayRange',max(max(converted_to_half{:,:})),'Symmetric',false,...
%     'RowLabels',converted_to_half.Properties.RowNames,'ColumnLabels',...
%     converted_to_half.Properties.RowNames);
% 
% ztest2 = clustergram(pct_raw_table{:,:},'Colormap',parula_plus,...
%     'DisplayRange',max(max(pct_raw_table{:,:})),'Symmetric',false,...
%     'RowLabels',pct_raw_table.Properties.RowNames,'ColumnLabels',...
%     pct_raw_table.Properties.RowNames); 
% 
% make_connectivity_heatmap(pct_raw_table,false);
% xtickangle(90)
% title("as percent")
% colormap(parula_plus)
%% now approach similiary to before
% add a moa column
bidirectional_table.NEIGHBOR_MOA = change_to_moa(bidirectional_table.NEIGHBOR,all_categories);

%% from earlier to make the old moa table, maybe use it?

% create empty moa table 
bi_moa_table  = cell2table(cell(0,3), 'VariableNames', {'DRUG', 'NEIGHBOR', 'CONNECTIONS'});


% loop through each drug 
for i = drugs
    drug = i{1};
    
    % get just rows for that drug
    row_inds = find(strcmp(bidirectional_table.DRUG,drug));
    
    drug_rows = bidirectional_table(row_inds,:);
    
    % now want to go through each unique MoA, weighting with the weights
    % that are weight
    moas = unique(drug_rows.NEIGHBOR_MOA)';
    
    for j = moas
        moa = j{1};
        
        moa_row_inds = find(strcmp(drug_rows.NEIGHBOR_MOA,moa));
        moa_rows = drug_rows(moa_row_inds,:);
        
        % get sum of weights
        sum_of_weights = sum(moa_rows.CONNECTIONS);
        
        % make a new row for the table bby
        new_row = table({drug},{moa},sum_of_weights, 'VariableNames', {'DRUG', 'NEIGHBOR', 'CONNECTIONS'});
        
        bi_moa_table =vertcat(bi_moa_table,new_row);
    end 
    
end

%% want to make a percent one in this format for ease of readability



%% now we've made bidrectional moa table, can create that heatmap 
% though still need to scale by number of connections for each 

% resort by moa
bi_moa_table = sortrows(bi_moa_table,'NEIGHBOR');

all_moas = unique(bi_moa_table.NEIGHBOR)';
% allows grouping by MoA
all_drugs = unique(results_table.DRUG,'stable')';

empty_array = zeros(length(unique(bi_moa_table.NEIGHBOR)),length(unique(bi_moa_table.DRUG)));

two_way_table_bi_moa = array2table(empty_array,'VariableNames',unique(results_table.DRUG,'stable'),...
    'RowNames',unique(bi_moa_table.NEIGHBOR));

for i = all_drugs 
    drug = i{1};
    
    row_inds = find(strcmp(bi_moa_table.DRUG,drug));
    
    drug_rows = bi_moa_table(row_inds,:);
    
    neighbor_list = drug_rows.NEIGHBOR;
    for j = all_moas
        moa = j{1};
        if ismember(moa,neighbor_list)
            % find where it is in the drug rows
            n_ind = find(strcmp(drug_rows.NEIGHBOR,moa));
            weight = drug_rows.CONNECTIONS(n_ind);
            two_way_table_bi_moa{moa,drug} = weight;
        else
            two_way_table_bi_moa{moa,drug} = 0;
        end 
        
    end 
    
    
end 

%% now that we have this table, need to scale by num connections
% make another table so we can still have a table with the raw valuese 
two_way_table_bi_moa_pct = two_way_table_bi_moa;
% divide each by the total number 
two_way_table_bi_moa_pct{:,:} = two_way_table_bi_moa{:,:}./sum(two_way_table_bi_moa{:,:})*100;


%% generate heatmap 

if ~crossval_mode
    make_connectivity_heatmap(two_way_table_bi_moa_pct,true);
    colormap(parula_plus)
    title({bayes_file,"bidirectional"},'Interpreter','none')
end 
% show the number of connections for each in a table 
all_drugs = two_way_table_bi_moa_pct.Properties.VariableNames';
all_connection_sums = sum(two_way_table_bi_moa{:,:})';

connection_table = table(all_drugs,all_connection_sums,'VariableNames',{'DRUG','TOTAL_CONNECTIONS'});


%% want to loop through this moa table to see if things are correctly classified

all_drugs = two_way_table_bi_moa_pct.Properties.VariableNames;
all_moas =  two_way_table_bi_moa_pct.Properties.RowNames';

% want to loop through all the drugs then all the moas

max_moas = cell(length(all_drugs),1);

ct = 1;
for i = all_drugs
    drug = i{1};
    highest_val = 0;
    max_moa = {};
    for j = all_moas
        moa = j{1};
        
        % percent for this moa
        pct_target = two_way_table_bi_moa_pct{moa,drug};
        
        % if it's the highest, keep it
        if pct_target > highest_val
            highest_val = pct_target;
            max_moa = moa;
        end 
        
        
    end 
    
    max_moas{ct} = max_moa;
    ct = ct+1;
end 

% make a new table with the max moas 

drug_moas = change_to_moa(all_drugs,all_categories)';


% get percent success by comparing drug_moas and max_moas
all_correct = 0;
for i = 1:length(drug_moas)
    
    if strcmp(drug_moas{i},max_moas{i})
        all_correct = all_correct + 1;
    end 
        
    
end 

success_pct = round(all_correct/length(max_moas)*100)

% put in table
the_table = table(drug_moas,max_moas,'VariableNames',{'DRUG_MOA','MOA_MATCH'},'RowNames',all_drugs);

%remove untreated and water and calculate the success w/o them
try 
    no_unt_and_water = the_table;

    no_unt_and_water('Untreated',:) =[];
    no_unt_and_water('water',:) = [];

    % get percent success by comparing drug_moas and max_moas
    all_correct_no_control = 0;
    for i = 1:length(no_unt_and_water.DRUG_MOA)

        if strcmp(no_unt_and_water.DRUG_MOA{i},no_unt_and_water.MOA_MATCH{i})
            all_correct_no_control = all_correct_no_control + 1;
        end 


    end 

    success_pct_no_control = round(all_correct_no_control/length(no_unt_and_water.DRUG_MOA)*100)

catch
    disp("no untreated and water to remove")
end 
%% full raw table but scaled by column
% like the moa one where each column adds up to 100% based on the number of
% connections each ended up happening 
bi_raw_col_scaled = full_raw_table;

bi_raw_col_scaled{:,:} = full_raw_table{:,:}./sum(full_raw_table{:,:})*100;

%% maybe search through two way table bi moa pct --> pick out highest value 
all_moas = two_way_table_bi_moa_pct.Properties.RowNames';
% also want to go through symmetric_table to get the nn drug 
%bi_raw_col_scaled is scaled at 100% for each drug--this makes the most
%sense since it's the same method as the MoA 
all_drugs = all_drugs';
if crossval_mode
    for q = all_drugs'
        current_drug = q{1};

        % go through each moa and find highest value
        highest_moa_val = 0;
        match_moa = "";
        for j = all_moas
            current_moa = j{1};
            pct_val = two_way_table_bi_moa_pct{current_moa,current_drug};
            if pct_val > highest_moa_val
                highest_moa_val = pct_val; 
                match_moa = current_moa;
            end 

        end 
        
        % now loop through all drugs and use bi_raw_col_scaled table
        highest_drug_val = 0;
        match_drug = "";
        
        % use symemtric table since that's what we make the heatmap for,
        % but I feel like it makes less sense here since we aren't looking
        % as much at overarching and are more focused on the connections
        % that just that drug has
        highest_drug_val_symmetric = 0;
        match_drug_symmetric = "";
        for w = all_drugs'
           current_match_drug = w{1}; 
           pct_val = bi_raw_col_scaled{current_match_drug,current_drug};
           
           % probably won't use this value, but grabbing it just in case 
           pct_val_symmetric = symmetric_table{current_match_drug,current_drug};
           if pct_val > highest_drug_val
               highest_drug_val = pct_val;
               match_drug = current_match_drug;
           end 
            
           if pct_val_symmetric > highest_drug_val_symmetric
               highest_drug_val_symmetric = pct_val_symmetric;
               match_drug_symmetric = current_match_drug;
           end 
        end
        
       
        disp_str = strcat("Highest MoA for ", current_drug, " is ", match_moa, " at ",...
            num2str(round(highest_moa_val)), "%");
        disp_str2 = strcat("... and the highest drug is ", match_drug," at ", num2str(round(highest_drug_val)),"%");
        disp_str3 = strcat("... and the highest drug is ", match_drug_symmetric, ...
            " at " , num2str(round(highest_drug_val_symmetric)),"% for symmetric");
        disp(disp_str)
        disp(disp_str2)
        disp(disp_str3)
        disp("--")
    end 

end 