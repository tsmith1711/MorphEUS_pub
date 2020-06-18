%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MorpHEUS using MicrobeJ                                 %
% Created by Ian Richardson, Sophia Hu and Michaela Olson %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Description:
% This file takes raw MicrobeJ data files:
%   1. Separates bacteria based on which drug (folder) and image it came from
%   2. Filters out bacteria according to thresholds on its transverse profile
%   3. Filters out bacteria according to MicrobeJ documentation
%   4. Takes the average value of each variable for each image
%   5. Takes the CV and skewness of of each variable for each image
%   6. Creates a final data table and normalizes each variable 
%   7. Calculates principal components
%   8. Plots the 3D and 2D PCAs and Loadings
%   9. Optionally additonally plots 3D and 2D LDA plot and rugplots
%   10. Auto saves workspace

clear all 

%% Import libraries and fuctions

addpath(genpath('./segmentation_functions/'));
addpath(genpath('./lib/'));

% save name of git repo 
repo_folder = pwd;

% select parent directory with data and function subfolders 
dname = uigetdir();
cd(dname)

% add path to data sub folder
data_path = strcat(dname,'/data/');
addpath(genpath(data_path));

disp(strcat("working in directory ", dname))

%% Import Settings from YAML config
MINIMUM_REQUIRED_CONFIG_VERSION = 1.6;

% if user has a local config they want to use, do that
try
    config = ReadYaml(strcat(dname,'/config.yml'));
    disp("Using local config file")
catch
    % if not, use the default yaml file
    config = ReadYaml(strcat(repo_folder,'/config.yml'));
    disp("Using default config file")
end 

%ensure version of config file is greater than the minimum required verison
if config.metadata.version < MINIMUM_REQUIRED_CONFIG_VERSION
    error(strcat('Config version too old! Version is ', config.metadata.version, ' while minimum required version for this script has been specified as ', MINIMUM_REQUIRED_CONFIG_VERSION));
end

%% User settings

% create cell or image workspace
cell_or_img = 'Image';

% merge technical replicates
merge_tech_reps = true;

%% Does the user want to use the cell or the image workspace?

% if not hard coded, prompt user
if ~exist('cell_or_img','var')
    cell_or_img = questdlg('Image or cell workspace?', ...
        '?????', ...
        'Image','Cell','Cell');

    % if user clicks x, exit script 
    if isempty(cell_or_img)
        return
    end
end 

switch cell_or_img
case 'Cell'
    cell_mode = true;
case 'Image'
    cell_mode = false;
end
%% Channel selection -- skipped for now, we are using all channels 
channel_list = {'Channel 1 - Bacteria','Channel 2 - syto','Channel 3 - FM'};
%[indx,choice] = listdlg('ListString',channel_list,'InitialValue',[1 2 3]);

chosen_channels = channel_list;

% if user selected Channel 1, set bacteria to true
if any(strcmp(chosen_channels,'Channel 1 - Bacteria'))
    bacteria_tf = true;
end 

% if user selected Channel 2, set syto to true
if any(strcmp(chosen_channels,'Channel 2 - syto'))
    syto_tf = true;
end 

% if user selected Channel 3, set FM to true
if any(strcmp(chosen_channels,'Channel 3 - FM'))
    FM_tf = true;
end 


%% Start the timer
tic %start timer
fprintf('Timer started \n')

%% load .csv files

%supress meaningless warnings
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
warning('off', 'MATLAB:textscan:UnableToGuessFormat');

%selects specific columns from bacteria/fluorescence/transverse files to keep in the table 
bacteria_column_names = [{'NAME','NAME_id','IMAGE_meta','EXPERIMENT','EXPERIMENT_fullname'} config.bacteria_cols];
fluorescence_column_names = [{'NAME','NAME_id','IMAGE_meta', 'PARENT','PARENT_localization'} config.fluorescence_cols];
transverse_column_names = {'NAME','NAME_id','INDEX_rel','INTENSITY_ch1'};

%contains channel 1 bacteria cells with SHAPE attributes
bacteria_file = csv_combine(data_path,'bacteria',bacteria_column_names);

if FM_tf
    %contains channel 3 FM fluorescence data with SHAPE attributes 
    FM_fluorescence_file = csv_combine(data_path,'feature2',fluorescence_column_names);
end 

if syto_tf
    %contains channel 2 syto fluorescence data with SHAPE attributes 
    syto_fluorescence_file = csv_combine(data_path,'feature1',fluorescence_column_names);
end 

%contains profile of intensity values for each cell in dataset
transverse_file = csv_combine(data_path,'transverse',transverse_column_names);

% get numeric columns
numeric_bacteria_cols = varfun(@isnumeric,bacteria_file,'OutputFormat', 'uniform');

if FM_tf
    numeric_FM_cols = varfun(@isnumeric,FM_fluorescence_file,'OutputFormat', 'uniform');
end 
if syto_tf
    numeric_syto_cols = varfun(@isnumeric,syto_fluorescence_file,'OutputFormat', 'uniform');
end 

%turn warnings back on just in case
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
warning('on', 'MATLAB:textscan:UnableToGuessFormat');

fprintf('Finished uploading data. Total time elapsed: %.2fs \n',toc)

%% Remove photo channel header and file extension to allow match across files
%%%% 'c1/3 - Phase_BDQ_image_01_R3D.dv' -> 'Phase_BDQ_image_01_R3D'

%remove channel header
%remove file extension
bacteria_file.IMAGE_meta = cellfun(@(x)extractAfter(x, ' - '), bacteria_file.IMAGE_meta, 'UniformOutput', false);
for i = 1:size(bacteria_file.IMAGE_meta, 1)
    bacteria_file.IMAGE_meta{i}(find(bacteria_file.IMAGE_meta{i} == '.',1,'last'):end) = [];
end
if FM_tf
    FM_fluorescence_file.IMAGE_meta = cellfun(@(x)extractAfter(x, ' - '), FM_fluorescence_file.IMAGE_meta, 'UniformOutput', false);
    for i = 1:size(FM_fluorescence_file.IMAGE_meta, 1)
        FM_fluorescence_file.IMAGE_meta{i}(find(FM_fluorescence_file.IMAGE_meta{i} == '.',1,'last'):end) = [];
    end
end 
if syto_tf
    syto_fluorescence_file.IMAGE_meta = cellfun(@(x)extractAfter(x, ' - '), syto_fluorescence_file.IMAGE_meta, 'UniformOutput', false);
    for i = 1:size(syto_fluorescence_file.IMAGE_meta, 1)
        syto_fluorescence_file.IMAGE_meta{i}(find(syto_fluorescence_file.IMAGE_meta{i} == '.',1,'last'):end) = [];
    end
end 


fprintf('Simplified image filenames. Total time elapsed: %.2fs \n',toc)

%% Split bacteria and fluorescence tables by drug and photo
%%%%Possible bug~~photo may have nothing recognized in bacteria but may have fluorescence

%~~~~~~~~~~ find the drug and photo names ~~~~~~~~~~~~~~~%

%%% find all drugs (names of folders which were batched into MicrobeJ)
drugs = unique(bacteria_file.EXPERIMENT).';
ids = unique(bacteria_file.EXPERIMENT_fullname).';

%find the unique photo names for each drug for each file type and convert to struct for better labeling
photos.bacteria = cellfun(@(x)unique(bacteria_file{strcmp(bacteria_file.EXPERIMENT, x), 'IMAGE_meta'}).', drugs, 'UniformOutput', false);
photos.bacteria = cell2struct(photos.bacteria, drugs, 2);
if FM_tf
    photos.FM = cellfun(@(x)unique(bacteria_file{strcmp(bacteria_file.EXPERIMENT, x), 'IMAGE_meta'}).', drugs, 'UniformOutput', false);
    photos.FM = cell2struct(photos.FM, drugs, 2);
end 
if syto_tf
    photos.syto = cellfun(@(x)unique(bacteria_file{strcmp(bacteria_file.EXPERIMENT, x), 'IMAGE_meta'}).', drugs, 'UniformOutput', false);
    photos.syto = cell2struct(photos.syto, drugs, 2);
end 

%get full list of recognized photos whether some are missing in fluorescence (expected) or bacteria (unlikely)
%a set of photos that will include all of them even if one doesnt have bacteria data, but does have fluorescence data
%use this list of photos when matching is important, but must provide a check that the file actually has the image (check that the field has the struct)

if syto_tf & FM_tf
    for i = drugs
        drug = i{1};
        photos.master.(drug) = unique([photos.bacteria.(drug) photos.FM.(drug) photos.syto.(drug)]);
    end
elseif syto_tf
    for i = drugs
        drug = i{1};
        photos.master.(drug) = unique([photos.bacteria.(drug) photos.syto.(drug)]);
    end
elseif FM_tf
    for i = drugs
        drug = i{1};
        photos.master.(drug) = unique([photos.bacteria.(drug) photos.FM.(drug)]);
    end
else
    for i = drugs
        drug = i{1};
        photos.master.(drug) = unique([photos.bacteria.(drug)]);
    end
end 

%~~~~~~~~~~~~~~~ Split files ~~~~~~~~~~~~ %

%split main files into drugs and photos
bacteria = structfun( @(x) ... %for each drug
               cellfun(@(y) ... %for each bacteria photo name
                    bacteria_file(strcmp(bacteria_file.IMAGE_meta, char(y)), :), ... %find bacteria rows with that photo name
               x, 'UniformOutput', false), ...
           photos.bacteria, 'UniformOutput', false);
       
if syto_tf
    syto_fluorescence = structfun( @(x) ... %for each drug
                            cellfun(@(y) ... %for each syto photo name
                                syto_fluorescence_file(strcmp(syto_fluorescence_file.IMAGE_meta, char(y)), :), ... %find syto rows with that photo name
                            x, 'UniformOutput', false), ...
                        photos.master, 'UniformOutput', false);
end 

if FM_tf
    FM_fluorescence = structfun( @(x) ... %for each drug
                           cellfun(@(y) ... %for each FM photo name
                                FM_fluorescence_file(strcmp(FM_fluorescence_file.IMAGE_meta, char(y)), :), ... %find FM rows with that photo name
                           x, 'UniformOutput', false), ...
                      photos.master, 'UniformOutput', false);
end 
              
%convert to structs with each field being an image 
for i = drugs
    drug = i{1};
    bacteria.(drug) = cell2struct(bacteria.(drug), photos.bacteria.(drug), 2);
    if syto_tf
        syto_fluorescence.(drug) = cell2struct(syto_fluorescence.(drug), photos.master.(drug), 2);
    end
    if FM_tf
        FM_fluorescence.(drug) = cell2struct(FM_fluorescence.(drug), photos.master.(drug), 2);
    end 
end

%set bacteria, syto, and FM rownames to their NAME column
for i = drugs
    drug = i{1};
    %bacteria
    for j = photos.bacteria.(drug)
        photo = j{1};
        bacteria.(drug).(photo).Properties.RowNames = bacteria.(drug).(photo).NAME;
    end
    %FM
    if FM_tf
        for j = photos.master.(drug)
            photo = j{1};
          %  FM_fluorescence.(drug).(photo).Properties.RowNames = FM_fluorescence.(drug).(photo).NAME;
        end
    end 
    %syto
    if syto_tf
        for j = photos.master.(drug)
            photo = j{1};
         %   syto_fluorescence.(drug).(photo).Properties.RowNames = syto_fluorescence.(drug).(photo).NAME;
        end
    end 
end

fprintf('Begin splitting transverse based on bacteria file images and by bacteria. Total time elapsed: %.2fs \n',toc)

%% Split transverse
%split transverse based on bacteria file images and by bacteria 
%this is where 99% of the time goes sad :((((

fprintf('%d different drugs to process \n',length(fieldnames(bacteria)))

transverse = structfun( @(x) printstructfun(x,transverse_file), ... %for each drug in bacteria
              bacteria,'UniformOutput', false); 
         
 fprintf('Finished splitting transverse based on bacteria file images and by bacteria. Total time elapsed: %.2fs \n',toc)
         
%split bacteria into struct for transverse
for i = drugs
    drug = i{1};
    for j = photos.bacteria.(drug)
        photo = j{1};
        % Identifies and removes any missing bacteria profiles ( empty cells in transverse.(drug).(photo) )
        missing_profiles = find(cellfun(@isempty, transverse.(drug).(photo)));
        transverse.(drug).(photo)(missing_profiles) = [];
        
        %find list of bacteria per photo
        bact_list = cellfun(@(x) x.NAME{1}, transverse.(drug).(photo), 'UniformOutput', false);
        %convert to struct with bacteria names as fields
        transverse.(drug).(photo) = cell2struct(transverse.(drug).(photo), bact_list, 2);
    end
end

%% Filters out blurry + bad bacteria
%This section iterates through each bacterium profile in the transverse file.
%The transverse file contains a list of intensities across the transverse
%axis of a cell, and has been split into a tree of tables containing the
%raw data of each bacterium.
%Then the range, slope, and fit of the bacteria is calculated.
%Bacteria that fall within a specified set of parameters will be added to
%the bacteria profile cell array.
%Start: Transverse tree
%End: Two trees of tables per image: one tree of tables in which each row represents the
%blur data on a 'good' bacterium, and one tree of tables in which each row
%represents a 'bad' bacterium.

%%%%%%%%%%%% !!!IMPORTANT!!! vv if change also change row add at end of loop below
blur_col_names = {'NAME','range_1','range_2','range_avg','slope_1','slope_2','slope_avg','num_pks','sw_p_val', 'min_mid_index_diff', 'sw_p_val_match', 'maxima_match', 'range_match', 'range_diff_match', 'slope_match', 'slope_diff_match', 'min_mid_match'};
%%%%%%%%%%%%%%%%%%%%%%%%%

%structure for data on passed bacteria
blur_data = structfun(@(x) ...
                structfun(@(y) ... 
                    cell2table(cell(0,size(blur_col_names,2)), 'VariableNames', blur_col_names), ...
                x, 'UniformOutput', false), ...
            bacteria, 'UniformOutput', false);

%structure for data on failed bacteria
bad_bacteria_blur_data = structfun(@(x) ...
                structfun(@(y) ... 
                    cell2table(cell(0,size(blur_col_names,2)), 'VariableNames', blur_col_names), ...
                x, 'UniformOutput', false), ...
            bacteria, 'UniformOutput', false);

        
fprintf('Filtered out bad or blurry bacteria. Total time elapsed: %.2fs \n',toc)
%% Add any bacteria which don't appear in transverse to bad_bacteria_blur_data
for i = drugs
    drug = i{1};
    for j = photos.bacteria.(drug)
        photo = j{1};
        %find vector of missing profiles in image
        missing_bacteria = ~ismember( bacteria.(drug).(photo).NAME, fieldnames(transverse.(drug).(photo)) );
        
        %add as row(s) to bad_bacteria_blur_data table
        if sum(missing_bacteria) > 0
            %create new rows for each missing bacterium in the image
            rows_data = cell(sum(missing_bacteria), size(blur_col_names,2)-1);
            rows_data(:) = {'Missing_profile'};
            new_rows = [bacteria.(drug).(photo).NAME(missing_bacteria, :), rows_data];
            %add the rows to the table
            bad_bacteria_blur_data.(drug).(photo) = [bad_bacteria_blur_data.(drug).(photo) ; new_rows];
        end
   end
end

clear missing_bacteria rows_data new_rows;   
%% If any nan in a bacteria profile remove from transverse and add to bad_bacteria_blur_data
transverse_no_nan = transverse;

for i = drugs
    drug = i{1};
    for j = photos.bacteria.(drug)
        photo = j{1};
        for k = fieldnames(transverse.(drug).(photo)).'
            bacterium = k{1};
            temp = transverse.(drug).(photo).(bacterium){:,{'INDEX_rel', 'INTENSITY_ch1'}};
            if sum(isnan(temp(:))) ~= 0 % if there are any nan
                %add to bad bacteria
                bad_bacteria_blur_data.(drug).(photo) = [bad_bacteria_blur_data.(drug).(photo) ; [{bacterium} num2cell(NaN([1 (size(blur_col_names,2)-1)]))] ];
                %remove from transverse
                transverse_no_nan.(drug).(photo) = rmfield(transverse_no_nan.(drug).(photo), bacterium);
            end
        end
    end
end
clear i j k drug photo bacterium temp
%% Find blur data and add to tables

%iterates through each unique bacteria in transverse file
for i = drugs % i iterates through drugs
    drug = i{1};
    for j = photos.bacteria.(drug) %j iterates through photos
        photo = j{1};
        for k = fieldnames(transverse_no_nan.(drug).(photo)).' %k iterates through bacteria names
            bacterium = k{1};
            
            %%% DEBUG %%%
            % drug      %
            % photo     %
            % bacterium %
            %%%%%%%%%%%%%
            
            %store the current bacterium profile in a nicer name
            profile = transverse_no_nan.(drug).(photo).(bacterium);
            
            %find number of points in profile
            num_points = size(profile, 1);

            %DEBUG: plots the x and y intensities for that cell
           % plot(profile.INDEX_rel, profile.INTENSITY_ch1, 'bo');

            %find the absolute min and middle point of the intensity profile
            cell_min = min(profile.INTENSITY_ch1);
            cell_min_index = find(profile.INTENSITY_ch1 == cell_min,1);
            cell_mid_index = round(num_points/2);
            min_mid_index_diff = abs(cell_min_index-cell_mid_index);

            %finds the absolate max of each side of the min
            cell_max_1 = max(profile.INTENSITY_ch1(1:cell_min_index));
            cell_max_2 = max(profile.INTENSITY_ch1(cell_min_index:num_points));
            cell_max_1_index = find(profile.INTENSITY_ch1(1:cell_min_index) == cell_max_1, 1);
            cell_max_2_index = find(profile.INTENSITY_ch1(cell_min_index:num_points) == cell_max_2, 1) + cell_min_index - 1; %find index in second half of profile and account for shift in index

            %calculates the range between cell_max and cell_min of both halves
            range_1 = cell_max_1-cell_min;
            range_2 = cell_max_2-cell_min;
            range_avg =(range_1+range_2)/2;
            range_diff = abs(range_1-range_2);

            %calculate the 25% and 75% intensity values
            y_lq_1 = range_1*0.25 + cell_min;
            y_uq_1 = range_1*0.75 + cell_min;
            y_lq_2 = range_2*0.25 + cell_min;
            y_uq_2 = range_2*0.75 + cell_min;

            %interpolate the INDEX_rel x values for the 25% and 75% intensity values
            x_uq_1 = interpolate_x(y_uq_1, profile.INDEX_rel(cell_max_1_index:cell_min_index).', profile.INTENSITY_ch1(cell_max_1_index:cell_min_index).');
            x_lq_1 = interpolate_x(y_lq_1, profile.INDEX_rel(cell_max_1_index:cell_min_index).', profile.INTENSITY_ch1(cell_max_1_index:cell_min_index).');
            x_uq_2 = interpolate_x(y_uq_2, profile.INDEX_rel(cell_min_index:cell_max_2_index).', profile.INTENSITY_ch1(cell_min_index:cell_max_2_index).');
            x_lq_2 = interpolate_x(y_lq_2, profile.INDEX_rel(cell_min_index:cell_max_2_index).', profile.INTENSITY_ch1(cell_min_index:cell_max_2_index).');
            
            %calculate slope between uq_1(75%) and lq_1(25%) intensity values
            slope_1 = abs((y_uq_1-y_lq_1)/(x_uq_1-x_lq_1));
            slope_2 = abs((y_uq_2-y_lq_2)/(x_uq_2-x_lq_2));
            slope_avg = (slope_1+slope_2)/2;
            slope_diff = abs(slope_1-slope_2);
            
            %find local maxima between sides
            if abs(cell_max_1_index-cell_max_2_index) > 3
                pks = findpeaks(profile.INTENSITY_ch1(cell_max_1_index:cell_max_2_index).');
                num_pks = length(pks);
            else
                pks = -1; %if there are only three or less points beteween maxima, must be bad cell
            end

            %calculate goodness of fit to gaussian
            sw_p_val = -1; 
            adjusted_curve = abs(profile.INTENSITY_ch1(cell_max_1_index:cell_max_2_index)-max(profile.INTENSITY_ch1));
            if abs(cell_max_1_index-cell_max_2_index) > 4 %swtest only works if there are more than 4 points 
                warning('off', 'all');
                [~, sw_p_val] = swtest(adjusted_curve); %shapiro-Wilk normality test
                warning('on', 'all');
            end

            %see which variables meet the constraints    
            sw_p_val_match = sw_p_val > config.blur_thresholds.sw_p_val_limit;
            maxima_match = num_pks <= config.blur_thresholds.maxima_limit;
            range_match = range_1 >= config.blur_thresholds.range_min && range_2 >= config.blur_thresholds.range_min;
            range_diff_match = range_diff < config.blur_thresholds.range_diff_max;
            slope_match = slope_1 > config.blur_thresholds.slope_min && slope_2 > config.blur_thresholds.slope_min;
            slope_diff_match = slope_diff < config.blur_thresholds.slope_diff_max;
            min_mid_match = min_mid_index_diff <= config.blur_thresholds.min_mid_diff_max;

            %Pick attributes to use for filter
            filter = eval(config.blur_thresholds.attributes_to_check);


            %adds blur attributes to bacteria profile cell array (if it falls
            %within certain requirements), or to bad bacteria profile cell array if
            %it doesnt
            
            %%%%%%%%%%% IMPORTANT!!! vv if change also change blur_data init at beginning of loop
            new_row = {bacterium,range_1,range_2,range_avg,slope_1,slope_2,slope_avg,num_pks,sw_p_val, min_mid_index_diff, sw_p_val_match, maxima_match, range_match, range_diff_match, slope_match, slope_diff_match, min_mid_match}; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if filter
                %add as row to blur data table
                blur_data.(drug).(photo) = [blur_data.(drug).(photo) ; new_row];
            else    
                %adds bacteria that do not pass the blur test to the bad profiles
                %cell array and remove from bacteria array
                bad_bacteria_blur_data.(drug).(photo) = [bad_bacteria_blur_data.(drug).(photo) ; new_row];
            end
        end
    end
end

fprintf('Continued work with blur data. Total time elapsed: %.2fs \n',toc)
%clear unnecessesary loop variables to get rid of clutter
clear adjusted_curve bact_list cell_max_1 cell_max_1_index cell_max_2 cell_max_2_index cell_mid_index cell_min cell_min_index ...
      filter i j k maxima_match min_mid_index_diff min_mid_match new_row num_pks num_points pks profile range_1 range_2 range_avg ...
      range_diff range_diff_match range_match slope_1 slope_2 slope_avg slope_diff slope_diff_match slope_match sw_p_val sw_p_val_match ...
      x_lq_1 x_lq_2 x_uq_1 x_uq_2 y_lq_1 y_lq_2 y_uq_1 y_uq_2 drug photo bacterium blur_thresholds


%% Remove filtered bacteria from bacteria and fluorescence files
%%%%in: bacteria, FM_fluorescence, syto_fluorescence
%%%%out:filtered_bacteria, filtered_FM_fluorescence,
%%%%    filtered_syto_fluorescence


%identify blurry/bad bacteria names from tree
bad_bacteria_names = structfun(@(x) ...
                         structfun(@(y) ...
                            y.NAME, ...
                         x, 'UniformOutput', false), ...
                     bad_bacteria_blur_data, 'UniformOutput', false);

%~~~~~~~~~~~~~~~~~~~ Remove from bacteria file ~~~~~~~~~~~~~~~~~%

%declare new file as direct copy of bacteria
filtered_bacteria = bacteria;

%remove all rows of bad bacteria
for i = drugs
    drug = i{1};
    for j = photos.bacteria.(drug)
        photo = j{1};
        if ~isempty(bad_bacteria_names.(drug).(photo)) %make sure there are bacteria to remove
            filtered_bacteria.(drug).(photo)(bad_bacteria_names.(drug).(photo), :) = []; 
        end
    end
end

%~~~~~~~~~~~~~~~~~~~ Remove from fluorescence files ~~~~~~~~~~~~~~~~~%

%%%%%%%%%% FM FLUORESCENCE %%%%%%
if FM_tf
    filtered_FM_fluorescence = FM_fluorescence;

    for i = drugs
        drug = i{1};
        for j = photos.master.(drug)
            photo = j{1};
            if ~isempty(bad_bacteria_names.(drug).(photo)) %make sure there are bacteria to remove
                %find binary matrix of rows whose parent is in bad_bacteria_names [1 = found(blurry), 0=not found(good)]
                badrows = cellfun(@(x) ismember(x, bad_bacteria_names.(drug).(photo)), FM_fluorescence.(drug).(photo).PARENT, 'UniformOutput', true);
                %remove all rows which whose parents were filtered out
                filtered_FM_fluorescence.(drug).(photo)(badrows, :) = [];
            end
        end
    end
    clear badrows
end 
%%%%%%%%%% SYTO FLUORESCENCE %%%%%%%%
if syto_tf
    filtered_syto_fluorescence = syto_fluorescence;

    for i = drugs
        drug = i{1};
        for j = photos.master.(drug)
            photo = j{1};
            if ~isempty(bad_bacteria_names.(drug).(photo)) %make sure there are bacteria to remove
                %find binary matrix of rows whose parent is in bad_bacteria_names [1 = found(blurry), 0=not found(good)]
                badrows = cellfun(@(x) ismember(x, bad_bacteria_names.(drug).(photo)), syto_fluorescence.(drug).(photo).PARENT, 'UniformOutput', true);
                %remove all rows which whose parents were filtered out
                filtered_syto_fluorescence.(drug).(photo)(badrows, :) = [];
            end
        end
    end
    clear badrows
end 
%% ADD EXTRA VARIABLES
% percent bacteria with fluorescence? (add onto averaged row of photo the number of rows != 0 in maxima count divided by numrows in the original split file)
% percent of each type of localization in image? (split from bacteria file and make own file (remove localization from fluorescence file)
% percent cell covered in fluorescence (look for children of bacterium add their area and divide by area of bacterium, if no children 0)

%% Merge technical replicate images if they exist 

% if user prompted to merge technical replicates 
if merge_tech_reps
    % use external function merge_reps to merge technical replicatese
    filtered_bacteria = merge_reps(filtered_bacteria);
    filtered_FM_fluorescence = merge_reps(filtered_FM_fluorescence);
    filtered_syto_fluorescence = merge_reps(filtered_syto_fluorescence);
end 

%% Merge tech rep spots together
% based on information in pad_info.csv, merge any technical replicates that
% are on different pads

% confirm that pad_info.csv exists -- if not, skip this 
if isfile('pad_info.csv')
    info_table = readtable('pad_info.csv','ReadVariableNames',true,'ReadRowNames',true);

    % use external function spot_combine to combine according to information in
    % pad_info.csv 
    filtered_bacteria = spot_combine(filtered_bacteria,info_table);
    filtered_FM_fluorescence = spot_combine(filtered_FM_fluorescence,info_table);
    filtered_syto_fluorescence = spot_combine(filtered_syto_fluorescence,info_table);
    % allow to use xlsx file as well just to be nice
elseif isfile('pad_info.xlsx')
    info_table = readtable('pad_info.xlsx','ReadVariableNames',true,'ReadRowNames',true);

    % use external function spot_combine to combine according to information in
    % pad_info.csv 
    filtered_bacteria = spot_combine(filtered_bacteria,info_table);
    filtered_FM_fluorescence = spot_combine(filtered_FM_fluorescence,info_table);
    filtered_syto_fluorescence = spot_combine(filtered_syto_fluorescence,info_table);
    
end 

%% Group photos together
%%%% merge a few photos together so that there is more data per point in
%%%% the final pca.
%%%% mergeImages.m removes rownames!!
%%%% in: bacteria, FM_fluorescence, syto_fluorescence
%%%% out: grouped_bacteria, grouped_FM_fluorescence,grouped_syto_fluorescence

grouped_bacteria = cell2struct(cell(1,size(drugs,2)), drugs, 2);
if FM_tf
    grouped_FM_fluorescence = cell2struct(cell(1,size(drugs,2)), drugs, 2);
end 
if syto_tf
    grouped_syto_fluorescence = cell2struct(cell(1,size(drugs,2)), drugs, 2);
end

for i = drugs
    drug = i{1};
    
    %find number of photos in drug
    num_photos = length(fieldnames(filtered_bacteria.(drug)));
    
    %make sure group size doesn't exceed number of photos
    if config.group_size > num_photos
        error(strcat('Group size larger than number of photos in: ', drug, ' (', num_photos, ')'));
    end
    
    %find the photos that will go into each group
    num_extra_photos = mod(num_photos, config.group_size);
    round_num_photos = num_photos-num_extra_photos;
    round_indexes_matrix = reshape(1:round_num_photos,config.group_size,round_num_photos/config.group_size).';
    group_indexes = mat2cell(round_indexes_matrix,ones(1,round_num_photos/config.group_size), config.group_size).';
    if num_extra_photos ~= 0
        group_indexes{length(group_indexes)} = [group_indexes{length(group_indexes)} num_photos-num_extra_photos+1:num_photos];
    end
    
    %use external function mergeImages to group photos 
    fnames = fieldnames(filtered_bacteria.(drug));
    grouped_bacteria.(drug) = cellfun(@(indexes) mergeImages(filtered_bacteria.(drug),fnames(indexes)), group_indexes, 'UniformOutput', false);
    
    if FM_tf 
        fnames = fieldnames(filtered_FM_fluorescence.(drug));
        grouped_FM_fluorescence.(drug) = cellfun(@(indexes) mergeImages(filtered_FM_fluorescence.(drug),fnames(indexes)), group_indexes, 'UniformOutput', false);
    end 
    if syto_tf
        fnames = fieldnames(filtered_syto_fluorescence.(drug));
        grouped_syto_fluorescence.(drug) = cellfun(@(indexes) mergeImages(filtered_syto_fluorescence.(drug),fnames(indexes)), group_indexes, 'UniformOutput', false);
    end 
end
clear num_photos num_extra_photos round_num_photos round_indexes_matrix group_indexes

%% Make sure there are no NaNs / missing values <<unimplemented>>
%in: grouped_(bacteria, FM, syto)
%out: grouped_no_nan_(bacteria, FM, syto)

%%%%%%EXTRA IDEAS FOR REPLACING NaN%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%replace weird missing indicators with NaN
    %standardizeMissing(table, indicator);
    %%replace NaN with 0
    %grouped_no_nan_bacteria.(drug){j} = fillmissing(table, 'constant', 0); or nanmean() of column
    %%remove rows with NaN
    %grouped_no_nan_bacteria.(drug){j} = rmmissing(table);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%USES EXTERNAL FUNCTION replace_nan

%bacteria
grouped_no_nan_bacteria = structfun(@(x) ...
                               cellfun(@(y) ...
                                   replace_nan(y, numeric_bacteria_cols), ...
                               x, 'UniformOutput', false), ...
                          grouped_bacteria, 'UniformOutput', false);
%syto
if syto_tf
    grouped_no_nan_syto_fluorescence = structfun(@(x) ...
                                           cellfun(@(y) ...
                                               replace_nan(y, numeric_syto_cols), ...
                                           x, 'UniformOutput', false), ...
                                       grouped_syto_fluorescence, 'UniformOutput', false);
end 
%FM
if FM_tf
    grouped_no_nan_FM_fluorescence = structfun(@(x) ...
                                           cellfun(@(y) ...
                                               replace_nan(y, numeric_FM_cols), ...
                                           x, 'UniformOutput', false), ...
                                       grouped_FM_fluorescence, 'UniformOutput', false);
end 
%% Beginning to filter out data based on MicrobeJ documentation
% creating file trees to store the bad data 


%~~~~~~make empty trees for each file type~~~~%

bad_bac = structfun(@(x) ... %loops through drugs
                cellfun(@(y) ... %loops through groups
                    cell2table(cell(0,size(bacteria_column_names,2)), 'VariableNames', bacteria_column_names), ...
                x, 'UniformOutput', false), ...
            grouped_no_nan_bacteria, 'UniformOutput', false);
if FM_tf         
    bad_FM = structfun(@(x) ... %loops through drugs
                    cellfun(@(y) ... %loops through groups
                        cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                    x, 'UniformOutput', false), ...
                grouped_no_nan_FM_fluorescence, 'UniformOutput', false);
end 
if syto_tf
    bad_syto = structfun(@(x) ... %loops through drugs
                    cellfun(@(y) ... %loops through groups
                        cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                    x, 'UniformOutput', false), ...
                grouped_no_nan_syto_fluorescence, 'UniformOutput', false);
end                         
                               
%% Filter out bacteria with circ too high 

% - bacteria - %

for i = drugs % for each drug 
    drug = i{1};
    bad_ct = 1;
    for l = 1:length(grouped_no_nan_bacteria.(drug)) 
        I = find(grouped_no_nan_bacteria.(drug){l}.SHAPE_circularity > 1);  
        if ~isempty(I)
             % if bad values have been found.....
            bad_bac.(drug){bad_ct} = grouped_no_nan_bacteria.(drug){l}(I,:);
            bad_ct = bad_ct+1;
            % need to additionally remove any children 
            parents = grouped_no_nan_bacteria.(drug){l}(I,:).NAME;
            for k = 1:length(parents)
                child_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,parents{k});
                if any(child_match)
                    grouped_no_nan_FM_fluorescence.(drug){l}(child_match,:) = [];
                end 
            end
            for k = 1:length(parents)
                child_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,parents{k});
                if any(child_match)
                    grouped_no_nan_syto_fluorescence.(drug){l}(child_match,:) = [];
                end 
            end
            grouped_no_nan_bacteria.(drug){l}(I,:) = [];       
        end
    end
end


% - syto - % 
if syto_tf
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_syto_fluorescence.(drug)) 
            if size(grouped_no_nan_syto_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_syto_fluorescence.(drug){l}.SHAPE_circularity > 1);
                if ~isempty(I)
                     % if bad values have been found.....
                    bad_syto.(drug){bad_ct} = grouped_no_nan_syto_fluorescence.(drug){l}(I,:);
                    bad_ct = bad_ct+1;
                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_syto_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_FM_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end

                    grouped_no_nan_syto_fluorescence.(drug){l}(I,:) = [];
                end
            end 
        end
    end
end 

% - FM - % 
if FM_tf
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_FM_fluorescence.(drug)) 
            % find indicies of rows where circularity is greater than 1
            if size(grouped_no_nan_FM_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_FM_fluorescence.(drug){l}.SHAPE_circularity > 1);
                if ~isempty(I)
                     % if bad values have been found.....
                    % remove where circularity is greater than 1 
                    bad_FM.(drug){bad_ct} = grouped_no_nan_FM_fluorescence.(drug){l}(I,:);
                    bad_ct = bad_ct+1;
                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_FM_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_syto_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end
                    grouped_no_nan_FM_fluorescence.(drug){l}(I,:) = [];
                end
            end 
        end
    end
end 

%% Filter out bacteria with solidity too high               
% per MicrobeJ documentation, solidity cannot be greater than 1

% - bacteria - %
for i = drugs % for each drug 
    drug = i{1};
    bad_ct = 1;
    for l = 1:length(grouped_no_nan_bacteria.(drug)) 
        if size(grouped_no_nan_bacteria.(drug){l},1)
            I = find(grouped_no_nan_bacteria.(drug){l}.SHAPE_solidity > 1);  
            if ~isempty(I)
                 % if bad values have been found.....
                if isempty(bad_bac.(drug){bad_ct})
                    %if this table is empty, add this value
                    bad_bac.(drug){bad_ct} = grouped_no_nan_bacteria.(drug){l}(I,:);
                else 
                    %if the table already has values, concatenate the new values
                    bad_bac.(drug){bad_ct} = [bad_bac.(drug){bad_ct}; grouped_no_nan_bacteria.(drug){l}(I,:)];
                end
                bad_ct = bad_ct+1;
                % need to additionally remove any children 
                parents = grouped_no_nan_bacteria.(drug){l}(I,:).NAME;
                for k = 1:length(parents)
                    child_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,parents{k});
                    if any(child_match)
                        grouped_no_nan_FM_fluorescence.(drug){l}(child_match,:) = [];
                    end 
                end
                for k = 1:length(parents)
                    child_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,parents{k});
                    if any(child_match)
                        grouped_no_nan_syto_fluorescence.(drug){l}(child_match,:) = [];
                    end 
                end
                grouped_no_nan_bacteria.(drug){l}(I,:) = [];
            end
        end 
    end
end


% - syto - % 
if syto_tf 
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_syto_fluorescence.(drug)) 
            if size(grouped_no_nan_syto_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_syto_fluorescence.(drug){l}.SHAPE_solidity > 1);
                if ~isempty(I)
                    % if bad values have been found.....
                    if isempty(bad_syto.(drug){bad_ct})
                        %if this table is empty, add this value
                        bad_syto.(drug){bad_ct} = grouped_no_nan_syto_fluorescence.(drug){l}(I,:);
                    else 
                        %if the table already has values, concatenate the new values
                        bad_syto.(drug){bad_ct} = [bad_syto.(drug){bad_ct}; grouped_no_nan_syto_fluorescence.(drug){l}(I,:)];
                    end
                    bad_ct = bad_ct+1;
                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_syto_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_FM_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end
                    grouped_no_nan_syto_fluorescence.(drug){l}(I,:) = [];
                end 
            end
        end
    end
end 

% - FM - % 
if FM_tf
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_FM_fluorescence.(drug)) 
            % find indicies of rows where solidity is greater than 1
            if size(grouped_no_nan_FM_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_FM_fluorescence.(drug){l}.SHAPE_solidity > 1);
                if ~isempty(I)
                    % if there are values with solidity > 1 
                    % remove where solidity is greater than 1 
                    if isempty(bad_FM.(drug){bad_ct})
                        %if this table is empty, add this value
                        bad_FM.(drug){bad_ct} = grouped_no_nan_FM_fluorescence.(drug){l}(I,:);
                    else 
                        %if the table already has values, concatenate the new values 
                        bad_FM.(drug){bad_ct} = [bad_FM.(drug){bad_ct}; grouped_no_nan_FM_fluorescence.(drug){l}(I,:)];
                    end
                    bad_ct = bad_ct+1;

                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_FM_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_syto_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end
                    grouped_no_nan_FM_fluorescence.(drug){l}(I,:) = [];
                end
            end 
        end
    end
end 

%% Filter out bacteria with an intensity of 0

% - syto - % 
if syto_tf 
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_syto_fluorescence.(drug)) 
            if size(grouped_no_nan_syto_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_syto_fluorescence.(drug){l}.INTENSITY == 0);
                if ~isempty(I)
                    % if bad values have been found.....
                    if isempty(bad_syto.(drug){bad_ct})
                        %if this table is empty, add this value
                        bad_syto.(drug){bad_ct} = grouped_no_nan_syto_fluorescence.(drug){l}(I,:);
                    else 
                        %if the table already has values, concatenate the new values
                        bad_syto.(drug){bad_ct} = [bad_syto.(drug){bad_ct}; grouped_no_nan_syto_fluorescence.(drug){l}(I,:)];
                    end
                    bad_ct = bad_ct+1;
                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_syto_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_FM_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end
                    grouped_no_nan_syto_fluorescence.(drug){l}(I,:) = [];
                end 
            end
        end
    end
end 

% - FM - % 
if FM_tf
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_FM_fluorescence.(drug)) 
            % find indicies of rows where solidity is greater than 1
            if size(grouped_no_nan_FM_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_FM_fluorescence.(drug){l}.INTENSITY == 0);
                if ~isempty(I)
                    % if there are values with solidity > 1 
                    % remove where solidity is greater than 1 
                    if isempty(bad_FM.(drug){bad_ct})
                        %if this table is empty, add this value
                        bad_FM.(drug){bad_ct} = grouped_no_nan_FM_fluorescence.(drug){l}(I,:);
                    else 
                        %if the table already has values, concatenate the new values 
                        bad_FM.(drug){bad_ct} = [bad_FM.(drug){bad_ct}; grouped_no_nan_FM_fluorescence.(drug){l}(I,:)];
                    end
                    bad_ct = bad_ct+1;

                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_FM_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_syto_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end
                    grouped_no_nan_FM_fluorescence.(drug){l}(I,:) = [];
                end
            end 
        end
    end
end 
%% Filter out bacteria with aspectRatio < 1
% per MicrobeJ documentation, results should have an aspectRatio >= 1

% - bacteria - %
for i = drugs % for each drug 
    drug = i{1};
    bad_ct = 1;
    for l = 1:length(grouped_no_nan_bacteria.(drug)) 
        I = find(grouped_no_nan_bacteria.(drug){l}.SHAPE_aspectRatio < 1);  
        if ~isempty(I)
            % if bad values have been found.....
            if isempty(bad_bac.(drug){bad_ct})
                %if this table is empty, add this value
                bad_bac.(drug){bad_ct} = grouped_no_nan_bacteria.(drug){l}(I,:);
            else 
                %if the table already has values, concatenate the new values
                bad_bac.(drug){bad_ct} = [bad_bac.(drug){bad_ct}; grouped_no_nan_bacteria.(drug){l}(I,:)];
            end
            bad_ct = bad_ct+1;
            % need to additionally remove any children 
            parents = grouped_no_nan_bacteria.(drug){l}(I,:).NAME;
            for k = 1:length(parents)
                child_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,parents{k});
                if any(child_match)
                    grouped_no_nan_FM_fluorescence.(drug){l}(child_match,:) = [];
                end 
            end
            for k = 1:length(parents)
                child_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,parents{k});
                if any(child_match)
                    grouped_no_nan_syto_fluorescence.(drug){l}(child_match,:) = [];
                end 
            end
            grouped_no_nan_bacteria.(drug){l}(I,:) = [];
        end
    end
end


%- syto - % 
if syto_tf
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_syto_fluorescence.(drug)) 
            if size(grouped_no_nan_syto_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_syto_fluorescence.(drug){l}.SHAPE_aspectRatio < 1);
                if ~isempty(I)
                    % if bad values have been found.....
                    if isempty(bad_syto.(drug){bad_ct})
                        %if this table is empty, add this value
                        bad_syto.(drug){bad_ct} = grouped_no_nan_syto_fluorescence.(drug){l}(I,:);
                    else 
                        %if the table already has values, concatenate the new values
                        bad_syto.(drug){bad_ct} = [bad_syto.(drug){bad_ct}; grouped_no_nan_syto_fluorescence.(drug){l}(I,:)];
                    end
                    bad_ct = bad_ct+1;


                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_syto_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_FM_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end
                    grouped_no_nan_syto_fluorescence.(drug){l}(I,:) = [];
                end 
            end
        end
    end
end 


%- FM - %  
if FM_tf
    for i = drugs % for each drug 
        drug = i{1};
        bad_ct = 1;
        for l = 1:length(grouped_no_nan_FM_fluorescence.(drug)) 
            % find indicies of rows where AR is less than 1
            if size(grouped_no_nan_FM_fluorescence.(drug){l},1)
                I = find(grouped_no_nan_FM_fluorescence.(drug){l}.SHAPE_aspectRatio < 1);
                if ~isempty(I)
                    % remove bacteria where AR is less than 1 
                    if isempty(bad_FM.(drug){bad_ct})
                        %if this table is empty, add this value
                        bad_FM.(drug){bad_ct} = grouped_no_nan_FM_fluorescence.(drug){l}(I,:);
                    else 
                        %if the table already has values, concatenate the new values 
                        bad_FM.(drug){bad_ct} = [bad_FM.(drug){bad_ct}; grouped_no_nan_FM_fluorescence.(drug){l}(I,:)];
                    end
                    bad_ct = bad_ct+1;


                    % need to additionally remove bacteria parents 
                    children = grouped_no_nan_FM_fluorescence.(drug){l}(I,:).PARENT;
                    for k = 1:length(children)
                        parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                        if any(parent_match)
                            grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                        end 
                    end 
                    % need to additionally remove siblings 
                    for k = 1:length(children)
                        sibling_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,children{k});
                        if any(sibling_match)
                            grouped_no_nan_syto_fluorescence.(drug){l}(sibling_match,:) = [];
                        end 
                    end
                    grouped_no_nan_FM_fluorescence.(drug){l}(I,:) = [];
                end 
            end
        end
    end
end 


%% Create outlier-included plots for observation later

% - bacteria - %

bacteria_hist_out = hist_file(drugs,numeric_bacteria_cols, bacteria_column_names,grouped_no_nan_bacteria);

% - syto - % 
if syto_tf
    try 
        syto_hist_out = hist_file(drugs,numeric_syto_cols, fluorescence_column_names,grouped_no_nan_syto_fluorescence);
    catch
        disp("unable to produce outlier syto histogram")
    end 
end 
% -- FM -- %
if FM_tf 
    try
        FM_hist_out = hist_file(drugs,numeric_FM_cols, fluorescence_column_names,grouped_no_nan_FM_fluorescence);
    catch
        disp("unable to produce outlier FM histogram")
    end 
end 

%% Getting rid of outliers 

bad_bac_out = [];
bad_FM_out = [];
bad_syto_out = [];

% sets the tolerance level. How many x the stdev to keep
numdev = 4;

% - bacteria - % 
bacteria_colnames_with_maxes = bacteria_column_names(numeric_bacteria_cols);
bacteria_colnames = bacteria_colnames_with_maxes(3:length(bacteria_colnames_with_maxes));
for i = drugs % for each drug 
    drug = i{1};
    for l = 1:length(grouped_no_nan_bacteria.(drug))        
        for j = 1:length(bacteria_colnames)
            variable_title = bacteria_colnames{j}; %get the name of the variable
            mean_val = mean(grouped_no_nan_bacteria.(drug){l}.(variable_title)); %find the mean val for the variable
            var_val = sqrt(var(grouped_no_nan_bacteria.(drug){l}.(variable_title))); % find the stdev for the variable 
            upperbound = var_val*numdev + mean_val; % find the upper bound for the kept data
            lowerbound = mean_val - var_val*numdev; %find the lower bound for the kept data
            if upperbound ~= lowerbound 
                I = find(grouped_no_nan_bacteria.(drug){l}.(variable_title) > upperbound); % if more than numdev deviations above the mean...
                if ~isempty(I)
                    % need to additionally remove any children 
                    parents = grouped_no_nan_bacteria.(drug){l}(I,:).NAME;
                    for k = 1:length(parents)
                        child_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,parents{k});
                        if any(child_match)
                            grouped_no_nan_FM_fluorescence.(drug){l}(child_match,:) = [];
                        end 
                    end
                    for k = 1:length(parents)
                        child_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,parents{k});
                        if any(child_match)
                            grouped_no_nan_syto_fluorescence.(drug){l}(child_match,:) = [];
                        end 
                    end
                    grouped_no_nan_bacteria.(drug){l}(I,:) = []; %delete 
                end
                I = find(grouped_no_nan_bacteria.(drug){l}.(variable_title) < lowerbound); % if more than numdev deviations below the mean...
                if ~isempty(I)
                    % need to additionally remove any children 
                    parents = grouped_no_nan_bacteria.(drug){l}(I,:).NAME;
                    for k = 1:length(parents)
                        child_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,parents{k});
                        if any(child_match)
                            grouped_no_nan_FM_fluorescence.(drug){l}(child_match,:) = [];
                        end 
                    end
                    for k = 1:length(parents)
                        child_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,parents{k});
                        if any(child_match)
                            grouped_no_nan_syto_fluorescence.(drug){l}(child_match,:) = [];
                        end 
                    end
                    grouped_no_nan_bacteria.(drug){l}(I,:) = []; % delete
                end
            end
        end
    end
end

% - syto - % 
if syto_tf
    syto_colnames = fluorescence_column_names(numeric_syto_cols);
    for i = drugs % for each drug 
        drug = i{1};
        for l = 1:length(grouped_no_nan_syto_fluorescence.(drug))      
            if size(grouped_no_nan_syto_fluorescence.(drug){l},1)
                for j = 1:length(syto_colnames)
                    variable_title = syto_colnames{j};
                    mean_val = mean(grouped_no_nan_syto_fluorescence.(drug){l}.(variable_title));
                    var_val = sqrt(var(grouped_no_nan_syto_fluorescence.(drug){l}.(variable_title)));
                    upperbound = var_val*numdev + mean_val;
                    lowerbound = mean_val - var_val*numdev;
                    I = find(grouped_no_nan_syto_fluorescence.(drug){l}.(variable_title) > upperbound);
                    if ~isempty(I)
                        % need to additionally remove bacteria parents 
                        children = grouped_no_nan_syto_fluorescence.(drug){l}(I,:).PARENT;
                        for k = 1:length(children)
                            parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                            if any(parent_match)
                                grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                            end 
                        end 
                        % need to additionally remove siblings 
                        for k = 1:length(children)
                            sibling_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,children{k});
                            if any(sibling_match)
                                grouped_no_nan_FM_fluorescence.(drug){l}(sibling_match,:) = [];
                            end 
                        end
                       grouped_no_nan_syto_fluorescence.(drug){l}(I,:) = [];


                    end
                    I = find(grouped_no_nan_syto_fluorescence.(drug){l}.(variable_title) < lowerbound);
                    if ~isempty(I)
                        % need to additionally remove bacteria parents 
                        children = grouped_no_nan_syto_fluorescence.(drug){l}(I,:).PARENT;
                        for k = 1:length(children)
                            parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                            if any(parent_match)
                                grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                            end 
                        end 
                        % need to additionally remove siblings 
                        for k = 1:length(children)
                            sibling_match = strcmp(grouped_no_nan_FM_fluorescence.(drug){l}.PARENT,children{k});
                            if any(sibling_match)
                                grouped_no_nan_FM_fluorescence.(drug){l}(sibling_match,:) = [];
                            end 
                        end
                        grouped_no_nan_syto_fluorescence.(drug){l}(I,:) = [];
                    end
                end
            end 
        end
    end
end 

% - FM - % 
if FM_tf
    FM_colnames = fluorescence_column_names(numeric_FM_cols);  
    for i = drugs % for each drug 
        drug = i{1};
        for l = 1:length(grouped_no_nan_FM_fluorescence.(drug)) 
            if size(grouped_no_nan_FM_fluorescence.(drug){l},1)
                for j = 1:length(FM_colnames)
                    variable_title = FM_colnames{j};
                    mean_val = mean(grouped_no_nan_FM_fluorescence.(drug){l}.(variable_title));
                    var_val = sqrt(var(grouped_no_nan_FM_fluorescence.(drug){l}.(variable_title)));
                    upperbound = var_val*numdev + mean_val;
                    lowerbound = mean_val - var_val*numdev;
                    I = find(grouped_no_nan_FM_fluorescence.(drug){l}.(variable_title) > upperbound);
                    if ~isempty(I)
                        % need to additionally remove bacteria parents 
                        children = grouped_no_nan_FM_fluorescence.(drug){l}(I,:).PARENT;
                        for k = 1:length(children)
                            parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                            if any(parent_match)
                                grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                            end 
                        end 
                        % need to additionally remove siblings 
                        for k = 1:length(children)
                            sibling_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,children{k});
                            if any(sibling_match)
                                grouped_no_nan_syto_fluorescence.(drug){l}(sibling_match,:) = [];
                            end 
                        end
                        grouped_no_nan_FM_fluorescence.(drug){l}(I,:) = [];
                    end
                    I = find(grouped_no_nan_FM_fluorescence.(drug){l}.(variable_title) < lowerbound);
                    if ~isempty(I)
                        % need to additionally remove bacteria parents 
                        children = grouped_no_nan_FM_fluorescence.(drug){l}(I,:).PARENT;
                        for k = 1:length(children)
                            parent_match = strcmp(grouped_no_nan_bacteria.(drug){l}.NAME,children{k});
                            if any(parent_match)
                                grouped_no_nan_bacteria.(drug){l}(parent_match,:) = [];
                            end 
                        end 
                        % need to additionally remove siblings 
                        for k = 1:length(children)
                            sibling_match = strcmp(grouped_no_nan_syto_fluorescence.(drug){l}.PARENT,children{k});
                            if any(sibling_match)
                                grouped_no_nan_syto_fluorescence.(drug){l}(sibling_match,:) = [];
                            end 
                        end
                        grouped_no_nan_FM_fluorescence.(drug){l}(I,:) = [];
                    end
                end
            end 
        end
    end
end 
%% Create files for histograms 
% uses helper function hist_file to create files ready to be transformed
% into histograms or rugplots
% see hist_file for further information

% - bacteria - %

bacteria_hist = hist_file(drugs,numeric_bacteria_cols, bacteria_column_names,grouped_no_nan_bacteria);

% - syto - % 
if syto_tf
    try 
        syto_hist = hist_file(drugs,numeric_syto_cols, fluorescence_column_names,grouped_no_nan_syto_fluorescence);
    catch
        disp("unable to produce syto histogram")
    end 
end 
% -- FM -- %
if FM_tf
    try 
        FM_hist = hist_file(drugs,numeric_FM_cols, fluorescence_column_names,grouped_no_nan_FM_fluorescence);
    catch
        disp("unable to produce FM histogram")
    end 
end 

%% Krista's code

% need to save grouped_no_nan files to get original skewness later
% make copy of those files to normal transform
if ~cell_mode 
    transformed_bacteria = grouped_no_nan_bacteria;

    transformed_FM = grouped_no_nan_FM_fluorescence;

    transformed_syto = grouped_no_nan_syto_fluorescence;
end 
    

%% For cell workspace, need to add parent zeros 
if cell_mode

    % - syto and FM - % 
    for i = drugs % for each drug 
        drug = i{1};
        for l = 1:length(grouped_no_nan_bacteria.(drug))        

            % get the parent bacteria from the FM and Syto and the list of all
            % possible parents from bacteria 
            parents = grouped_no_nan_bacteria.(drug){l}.NAME;
            if FM_tf
                parent_f = grouped_no_nan_FM_fluorescence.(drug){l}.PARENT;
                % get numeric and not numeric columns in FM
                numeric_f = varfun(@isnumeric,grouped_no_nan_FM_fluorescence.(drug){l},'OutputFormat', 'uniform');
                notnum_f = ~numeric_f;
                numeric_f_col_names = grouped_no_nan_FM_fluorescence.(drug){l}.Properties.VariableNames(numeric_f);
            end 
            if syto_tf
                parent_s = grouped_no_nan_syto_fluorescence.(drug){l}.PARENT;
                % get numeric and not numeric colimns in syto
                numeric_s = varfun(@isnumeric,grouped_no_nan_syto_fluorescence.(drug){l},'OutputFormat', 'uniform');
                notnum_s = ~numeric_s;
                numeric_s_col_names = grouped_no_nan_syto_fluorescence.(drug){l}.Properties.VariableNames(numeric_s);
            end 
            % get numeric and not numeric columns in bacteria
            numeric_b = varfun(@isnumeric,grouped_no_nan_bacteria.(drug){l},'OutputFormat', 'uniform');
            notnum_b = ~numeric_b;


            % check that there wasn't a bacteria removed earlier due to
            % blurriness that didn't also remove its child and remove the
            % corresponding child
            if syto_tf
                I = [];
                for z = 1:length(parent_s)
                    if ~any(strcmp(parent_s(z),parents))
                        I = [I z];
                    end
                end
                % removes syto children with no parent
                grouped_no_nan_syto_fluorescence.(drug){l}(I,:) = [];
            end 

            % removes FM children with no parent 
            if FM_tf
                I = [];
                for z = 1:length(parent_f)
                    if ~any(strcmp(parent_f(z),parents))
                        I = [I z];
                    end  
                end
                 grouped_no_nan_FM_fluorescence.(drug){l}(I,:) = [];
            end 
           % for all of the parent bacteria 
           for z = 1:length(parents)

               % find all parent bacteria that do not have a syto child
               if syto_tf
                    if ~any(strcmp(parents(z),parent_s))
                        % therefore, all the syto information should be 0
                        bac_info = grouped_no_nan_bacteria.(drug){l}(:,notnum_b);
                        bac_info2 = bac_info(z,:);
                        bac_info2 = bac_info2(:,1:3);
                        parent_info = grouped_no_nan_bacteria.(drug){l}.NAME;
                        parent_info2 = parent_info(z,:);
                        parent_info_row = table(parent_info2, "dummy", 'VariableNames',{'PARENT','PARENT_localization'});
                        empty_row = array2table(zeros(1,length(numeric_s_col_names)), 'VariableNames',numeric_s_col_names);

                        new_row = [bac_info2 parent_info_row empty_row];

                        % create a new syto row with all 0's for the given parent
                        grouped_no_nan_syto_fluorescence.(drug){l} = [grouped_no_nan_syto_fluorescence.(drug){l}; new_row];                

                    end
               end 
                % find all parent bacteria that do not have a FM child
                if FM_tf
                    if ~any(strcmp(parents(z),parent_f))
                        % therefore, all FM info should be zero
                        bac_info = grouped_no_nan_bacteria.(drug){l}(:,notnum_b);
                        bac_info2 = bac_info(z,:);
                        bac_info2 = bac_info2(:,1:3);
                        parent_info = grouped_no_nan_bacteria.(drug){l}.NAME;
                        parent_info2 = parent_info(z,:);
                        parent_info_row = table(parent_info2, "dummy", 'VariableNames',{'PARENT','PARENT_localization'});
                        empty_row = array2table(zeros(1,length(numeric_f_col_names)), 'VariableNames',numeric_f_col_names);
                        new_row = [bac_info2 parent_info_row empty_row];

                        % create a new FM row with all 0's for the given parent
                        grouped_no_nan_FM_fluorescence.(drug){l} = [grouped_no_nan_FM_fluorescence.(drug){l}; new_row];

                    end
                end 
            end

        end
    end

end 

%% median tables - only if in image workspace 
%%% get average for feature 1 and 2 though 
if ~cell_mode
    %~~~~~~make empty trees for each file type~~~~%

    center_bacteria_groups = structfun(@(x) ... %loops through drugs
                    cellfun(@(y) ... %loops through groups
                        cell2table(cell(0,size(bacteria_column_names,2)), 'VariableNames', bacteria_column_names), ...
                    x, 'UniformOutput', false), ...
                transformed_bacteria, 'UniformOutput', false);
    if FM_tf        
        center_FM_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_FM, 'UniformOutput', false);
    end 
    if syto_tf
        center_syto_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_syto, 'UniformOutput', false);
    end 

    %~~~~~~~~~~~~~BACTERIA~~~~~~~~~~~%
    numeric_cols = varfun(@isnumeric,bacteria_file,'OutputFormat', 'uniform'); %which cols contain numeric data
    feature_cols = contains(bacteria_file.Properties.VariableNames,'FEATURE');
    non_feature_numeric_cols = logical(numeric_cols - feature_cols);
    for i = drugs
        drug = i{1};
        for j = 1:size(transformed_bacteria.(drug), 2)
            %get table, drug name, and images of group
            atable = transformed_bacteria.(drug){j};
            experiment = atable(1, 'EXPERIMENT');
            experiment.Properties.RowNames = {};
            images = unique(atable.IMAGE_meta).';
            idname = unique(atable.EXPERIMENT_fullname);
            % just take first one
            idname = idname{1};
            %format images into nice table
            image_cols = cell(1, length(images));
            for k = 1:length(images)
                image_cols{k} = sprintf('IMAGE_%d',k);
            end
            photo_labels = cell2table(images, 'VariableNames', image_cols);
            exp_labels = cell2table({idname}, 'VariableNames', {'DRUG_id'});
            %horizontally merge experiment and image labels
            labels = [experiment exp_labels photo_labels];

            %average numeric data and add labels
            data = atable(:, non_feature_numeric_cols);
            avg_data = varfun(@nanmedian, data);
            %get averages for the feature data
            data = atable(:,feature_cols);
            avg_data2 = varfun(@nanmean, data);
            avg_data2.Properties.VariableNames = bacteria_column_names(feature_cols);
            avg_data.Properties.VariableNames = bacteria_column_names(non_feature_numeric_cols);
            
            center_bacteria_groups.(drug){j} = [labels avg_data2 avg_data];
        end
    end

    %~~~~~~~~~~~~~~~~FM~~~~~~~~~~~~~~~~~~~%
    if FM_tf
        numeric_cols = varfun(@isnumeric,FM_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        cols = fluorescence_column_names; %column names for file
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_FM_cols))), 'VariableNames', cols(numeric_FM_cols)); %numeric columns with one row of zero
        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_FM.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_FM.(drug){j};
                if size(atable)
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %average numeric data and add labels
                    data = atable(:, numeric_cols);
                    avg_data = varfun(@nanmedian, data);
                    avg_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    center_FM_groups.(drug){j} = [labels avg_data];        

                else
                    labels = cell(1, sum(~numeric_FM_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_FM_cols));
                    center_FM_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 
    %~~~~~~~~~~~~~SYTO~~~~~~~~~~~~~%
    if syto_tf
        cols = fluorescence_column_names; %column names for file
        numeric_cols = varfun(@isnumeric,syto_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_syto_cols))), 'VariableNames', cols(numeric_syto_cols)); %numeric columns with one row of zero

        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_syto.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_syto.(drug){j};
                if size(atable)
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %average numeric data and add labels
                    data = atable(:, numeric_cols);
                    avg_data = varfun(@nanmedian, data);
                    avg_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    center_syto_groups.(drug){j} = [labels avg_data];

                else
                    labels = cell(1, sum(~numeric_syto_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_syto_cols));
                    center_syto_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 

end 
clear image_cols images labels photo_labels table

%% Calculate Q1 

if ~cell_mode
    %~~~~~~make empty trees for each file type~~~~%

    q1_bacteria_groups = structfun(@(x) ... %loops through drugs
                    cellfun(@(y) ... %loops through groups
                        cell2table(cell(0,size(bacteria_column_names,2)), 'VariableNames', bacteria_column_names), ...
                    x, 'UniformOutput', false), ...
                transformed_bacteria, 'UniformOutput', false);
    if FM_tf        
        q1_FM_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_FM, 'UniformOutput', false);
    end 
    if syto_tf
        q1_syto_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_syto, 'UniformOutput', false);
    end 

    %~~~~~~~~~~~~~BACTERIA~~~~~~~~~~~%
    numeric_cols = varfun(@isnumeric,bacteria_file,'OutputFormat', 'uniform'); %which cols contain numeric data
    for i = drugs
        drug = i{1};
        for j = 1:size(transformed_bacteria.(drug), 2)
            %get table, drug name, and images of group
            atable = transformed_bacteria.(drug){j};
            experiment = atable(1, 'EXPERIMENT');
            experiment.Properties.RowNames = {};
            images = unique(atable.IMAGE_meta).';
            idname = unique(atable.EXPERIMENT_fullname);
            % just take first one
            idname = idname{1};
            %format images into nice table
            image_cols = cell(1, length(images));
            for k = 1:length(images)
                image_cols{k} = sprintf('IMAGE_%d',k);
            end
            photo_labels = cell2table(images, 'VariableNames', image_cols);
            exp_labels = cell2table({idname}, 'VariableNames', {'DRUG_id'});
            %horizontally merge experiment and image labels
            labels = [experiment exp_labels photo_labels];

            %average numeric data and add labels
            data = atable(:, numeric_cols);
            q1func = @(x) quantile(x,0.25);
            q1_data = varfun(q1func, data);
            q1_data.Properties.VariableNames = bacteria_column_names(numeric_cols);
            q1_bacteria_groups.(drug){j} = [labels q1_data];

        end
    end

    %~~~~~~~~~~~~~~~~FM~~~~~~~~~~~~~~~~~~~%
    if FM_tf
        numeric_cols = varfun(@isnumeric,FM_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        cols = fluorescence_column_names; %column names for file
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_FM_cols))), 'VariableNames', cols(numeric_FM_cols)); %numeric columns with one row of zero
        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_FM.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_FM.(drug){j};
                if size(atable)
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %average numeric data and add labels
                    data = atable(:, numeric_cols);
                    q1func = @(x) quantile(x,0.25);
                    q1_data = varfun(q1func, data);
                    q1_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    q1_FM_groups.(drug){j} = [labels q1_data];        
                else
                    labels = cell(1, sum(~numeric_FM_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_FM_cols));
                    q1_FM_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 
    %~~~~~~~~~~~~~SYTO~~~~~~~~~~~~~%
    if syto_tf
        cols = fluorescence_column_names; %column names for file
        numeric_cols = varfun(@isnumeric,syto_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_syto_cols))), 'VariableNames', cols(numeric_syto_cols)); %numeric columns with one row of zero

        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_syto.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_syto.(drug){j};
                if size(atable)
                  %  experiment = atable(1, 'EXPERIMENT');
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                  %  experiment.Properties.RowNames = {};
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %average numeric data and add labels
                    data = atable(:, numeric_cols);
                    q1func = @(x) quantile(x,0.25);
                    q1_data = varfun(q1func, data);
                    q1_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    q1_syto_groups.(drug){j} = [labels q1_data];

                else
                    labels = cell(1, sum(~numeric_syto_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_syto_cols));
                    q1_syto_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 

end 
clear image_cols images labels photo_labels table


%% Q3 

if ~cell_mode
    %~~~~~~make empty trees for each file type~~~~%

    q3_bacteria_groups = structfun(@(x) ... %loops through drugs
                    cellfun(@(y) ... %loops through groups
                        cell2table(cell(0,size(bacteria_column_names,2)), 'VariableNames', bacteria_column_names), ...
                    x, 'UniformOutput', false), ...
                transformed_bacteria, 'UniformOutput', false);
    if FM_tf        
        q3_FM_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_FM, 'UniformOutput', false);
    end 
    if syto_tf
        q3_syto_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_syto, 'UniformOutput', false);
    end 

    %~~~~~~~~~~~~~BACTERIA~~~~~~~~~~~%
    numeric_cols = varfun(@isnumeric,bacteria_file,'OutputFormat', 'uniform'); %which cols contain numeric data
    for i = drugs
        drug = i{1};
        for j = 1:size(transformed_bacteria.(drug), 2)
            %get table, drug name, and images of group
            atable = transformed_bacteria.(drug){j};
            experiment = atable(1, 'EXPERIMENT');
            experiment.Properties.RowNames = {};
            images = unique(atable.IMAGE_meta).';
            idname = unique(atable.EXPERIMENT_fullname);
            % just take first one
            idname = idname{1};
            %format images into nice table
            image_cols = cell(1, length(images));
            for k = 1:length(images)
                image_cols{k} = sprintf('IMAGE_%d',k);
            end
            photo_labels = cell2table(images, 'VariableNames', image_cols);
            exp_labels = cell2table({idname}, 'VariableNames', {'DRUG_id'});
            %horizontally merge experiment and image labels
            labels = [experiment exp_labels photo_labels];

            %average numeric data and add labels
            data = atable(:, numeric_cols);
            q3func = @(x) quantile(x,0.75);
            q3_data = varfun(q3func, data);
            q3_data.Properties.VariableNames = bacteria_column_names(numeric_cols);
            q3_bacteria_groups.(drug){j} = [labels q3_data];

        end
    end

    %~~~~~~~~~~~~~~~~FM~~~~~~~~~~~~~~~~~~~%
    if FM_tf
        numeric_cols = varfun(@isnumeric,FM_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        cols = fluorescence_column_names; %column names for file
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_FM_cols))), 'VariableNames', cols(numeric_FM_cols)); %numeric columns with one row of zero
        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_FM.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_FM.(drug){j};
                if size(atable)
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %average numeric data and add labels
                    data = atable(:, numeric_cols);
                    q3func = @(x) quantile(x,0.75);
                    q3_data = varfun(q3func, data);
                    q3_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    q3_FM_groups.(drug){j} = [labels q3_data];        
                else
                    labels = cell(1, sum(~numeric_FM_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_FM_cols));
                    q3_FM_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 
    %~~~~~~~~~~~~~SYTO~~~~~~~~~~~~~%
    if syto_tf
        cols = fluorescence_column_names; %column names for file
        numeric_cols = varfun(@isnumeric,syto_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_syto_cols))), 'VariableNames', cols(numeric_syto_cols)); %numeric columns with one row of zero

        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_syto.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_syto.(drug){j};
                if size(atable)
                  %  experiment = atable(1, 'EXPERIMENT');
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                  %  experiment.Properties.RowNames = {};
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %average numeric data and add labels
                    data = atable(:, numeric_cols);
                    q3func = @(x) quantile(x,0.75);
                    q3_data = varfun(q3func, data);
                    q3_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    q3_syto_groups.(drug){j} = [labels q3_data];

                else
                    labels = cell(1, sum(~numeric_syto_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_syto_cols));
                    q3_syto_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 

end 
clear image_cols images labels photo_labels table


%% IQR - only if in image workspace 

if ~cell_mode 
    %~~~~~~make empty trees for each file type~~~~%

    iqr_bacteria_groups = structfun(@(x) ... %loops through drugs
                    cellfun(@(y) ... %loops through groups
                        cell2table(cell(0,size(bacteria_column_names,2)), 'VariableNames', bacteria_column_names), ...
                    x, 'UniformOutput', false), ...
                transformed_bacteria, 'UniformOutput', false);

    if FM_tf        
        iqr_FM_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_FM, 'UniformOutput', false);
    end         

    if syto_tf
        iqr_syto_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_syto, 'UniformOutput', false);
    end 

    %~~~~~~~~~~~~~BACTERIA~~~~~~~~~~~%
    numeric_cols = varfun(@isnumeric,bacteria_file,'OutputFormat', 'uniform'); %which cols contain numeric data
    for i = drugs
        drug = i{1};
        for j = 1:size(transformed_bacteria.(drug), 2)
            %get table, drug name, and images of group
            atable = transformed_bacteria.(drug){j};
            experiment = atable(1, 'EXPERIMENT');
            experiment.Properties.RowNames = {};
            images = unique(atable.IMAGE_meta).';
            idname = unique(atable.EXPERIMENT_fullname);
            % just take first one
            idname = idname{1};
            %format images into nice table
            image_cols = cell(1, length(images));
            for k = 1:length(images)
                image_cols{k} = sprintf('IMAGE_%d',k);
            end
            photo_labels = cell2table(images, 'VariableNames', image_cols);
            %horizontally merge experiment and image labels
            labels = [experiment photo_labels];

            %get variance for numeric data and add labels
            data = atable(:, numeric_cols);
            var_data = varfun(@(x) iqr(x), data);
            var_data.Properties.VariableNames = bacteria_column_names(numeric_cols);
            iqr_bacteria_groups.(drug){j} = [labels var_data];
        end
    end        


    %~~~~~~~~~~~~~~~~FM~~~~~~~~~~~~~~~~~~~%
    if FM_tf
        numeric_cols = varfun(@isnumeric,FM_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        cols = fluorescence_column_names; %column names for file
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_FM_cols))), 'VariableNames', cols(numeric_FM_cols)); %numeric columns with one row of zero
        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_FM.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_FM.(drug){j};
                if size(atable)
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %get variance for numeric data and add labels
                    data = atable(:, numeric_cols);
                    var_data = varfun(@(x) iqr(x), data);
                    var_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    iqr_FM_groups.(drug){j} = [labels var_data];

                else
                    labels = cell(1, sum(~numeric_FM_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_FM_cols));
                    iqr_FM_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 

    %~~~~~~~~~~~~~SYTO~~~~~~~~~~~~~%
    if syto_tf
        cols = fluorescence_column_names; %column names for file
        numeric_cols = varfun(@isnumeric,syto_fluorescence_file,'OutputFormat', 'uniform'); %which cols contain numeric data
        zero_data = cell2table(num2cell(zeros(1,sum(numeric_syto_cols))), 'VariableNames', cols(numeric_syto_cols)); %numeric columns with one row of zero
        for i = drugs
            drug = i{1};
            for j = 1:size(transformed_syto.(drug), 2)
                %get table, drug name, and images of group
                atable = transformed_syto.(drug){j};
                if size(atable)
                    experiment = table({drug}, 'VariableNames',{'EXPERIMENT'});
                    images = unique(atable.IMAGE_meta).';

                    %format images into nice table
                    image_cols = cell(1, length(images));
                    for k = 1:length(images)
                        image_cols{k} = sprintf('IMAGE_%d',k);
                    end
                    photo_labels = cell2table(images, 'VariableNames', image_cols);

                    %horizontally merge experiment and image labels
                    labels = [experiment photo_labels];

                    %get variance for numeric data and add labels
                    data = atable(:, numeric_cols);
                    var_data = varfun(@(x) iqr(x), data);
                    var_data.Properties.VariableNames = fluorescence_column_names(numeric_cols);
                    iqr_syto_groups.(drug){j} = [labels var_data];
                else
                    labels = cell(1, sum(~numeric_syto_cols));
                    labels(:) = {'zeroData'};
                    labels = cell2table(labels, 'VariableNames', cols(~numeric_syto_cols));
                    iqr_syto_groups.(drug){j} = [labels zero_data];
                end 
            end
        end
    end 
end         

%% Adjust fluorescence column names and variance names 
%makes it easier to distinguish from regular bacteria column names
%adds a s_ in front of the syto column names and an f_ in front of fm column names

for i = drugs
    drug = i{1};
    % if FM was selected earlier 
    if FM_tf
        % if we are in cell mode 
        if cell_mode
            for j = 1:length(grouped_no_nan_FM_fluorescence.(drug))
                grouped_no_nan_FM_fluorescence.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('f_',x), grouped_no_nan_FM_fluorescence.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
            end
        else 
        % else we are in image mode 
            for j = 1:length(center_FM_groups.(drug))
                center_FM_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('f_',x), center_FM_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                %differentiate variance values 
                iqr_FM_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('f_',x,'_v'), iqr_FM_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                %... but make sure drug name is the same so we can join tables 
                iqr_FM_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(iqr_FM_groups.(drug){j}.Properties.VariableNames{1},'_v');
                %differentiate Q1 values
                q1_FM_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('f_',x,'_q1'), q1_FM_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                %... but make sure drug name is the same so we can join tables 
                q1_FM_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(q1_FM_groups.(drug){j}.Properties.VariableNames{1},'_q1');
                %differentiate Q3 values
                q3_FM_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('f_',x,'_q3'), q3_FM_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                %... but make sure drug name is the same so we can join tables 
                q3_FM_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(q3_FM_groups.(drug){j}.Properties.VariableNames{1},'_q3');
            end
        end
    end 
    % if syto was selected earlier 
    if syto_tf
        if cell_mode
            % if we are in cell mode
            for j = 1:length(grouped_no_nan_syto_fluorescence.(drug))
                grouped_no_nan_syto_fluorescence.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('s_',x), grouped_no_nan_syto_fluorescence.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
            end
        else 
            % else we are in image mode 
            for j = 1:length(center_syto_groups.(drug))
                center_syto_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('s_',x), center_syto_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                % differentiate variance values 
                iqr_syto_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('s_',x,'_v'), iqr_syto_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                % ... but make sure drug name is the same so we can join tables
                iqr_syto_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(iqr_syto_groups.(drug){j}.Properties.VariableNames{1},'_v');
                %differentiate Q1 values
                q1_syto_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('s_',x,'_q1'), q1_syto_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                %... but make sure drug name is the same so we can join tables 
                q1_syto_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(q1_syto_groups.(drug){j}.Properties.VariableNames{1},'_q1');
                %differentiate Q3 values
                q3_syto_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat('s_',x,'_q3'), q3_syto_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
                %... but make sure drug name is the same so we can join tables 
                q3_syto_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(q3_syto_groups.(drug){j}.Properties.VariableNames{1},'_q3');
            end
        end 
    end 
    if ~cell_mode
    % if not in cell mode, also need to adjust IQR and 
        for j = 1:length(center_bacteria_groups.(drug))
            iqr_bacteria_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat(x,'_v'), iqr_bacteria_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
             % ... but make sure drug name is the same so we can join tables
            iqr_bacteria_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(iqr_bacteria_groups.(drug){j}.Properties.VariableNames{1},'_v');
            
             q1_bacteria_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat(x,'_q1'), q1_bacteria_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
             % ... but make sure drug name is the same so we can join tables
             q1_bacteria_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(q1_bacteria_groups.(drug){j}.Properties.VariableNames{1},'_q1');
            
             q3_bacteria_groups.(drug){j}.Properties.VariableNames = cellfun(@(x) strcat(x,'_q3'), q3_bacteria_groups.(drug){j}.Properties.VariableNames, 'UniformOutput', false);
             % ... but make sure drug name is the same so we can join tables
            q3_bacteria_groups.(drug){j}.Properties.VariableNames{1} = extractBefore(q3_bacteria_groups.(drug){j}.Properties.VariableNames{1},'_q3');
        end
    end 
end


%% join average, variance and skew tables - only if img workspace

%~~~~~~make empty trees for each file type~~~~%
if ~cell_mode
    joined_bacteria_groups = structfun(@(x) ... %loops through drugs
                    cellfun(@(y) ... %loops through groups
                        cell2table(cell(0,size(bacteria_column_names,2)), 'VariableNames', bacteria_column_names), ...
                    x, 'UniformOutput', false), ...
                transformed_bacteria, 'UniformOutput', false);
    if FM_tf        
        joined_FM_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_FM, 'UniformOutput', false);
    end 
    if syto_tf
        joined_syto_groups = structfun(@(x) ... %loops through drugs
                        cellfun(@(y) ... %loops through groups
                            cell2table(cell(0,size(fluorescence_column_names,2)), 'VariableNames', fluorescence_column_names), ...
                        x, 'UniformOutput', false), ...
                    transformed_syto, 'UniformOutput', false);
    end         
    %~~~~~~~~~~~~~bacteria~~~~~~~~~~~~~%
    if ~cell_mode
        for i = drugs
            drug = i{1};
            for j = 1:size(iqr_bacteria_groups.(drug), 2)
                % join average table and variance table 
                joined_bacteria_groups.(drug){j} = join(center_bacteria_groups.(drug){j},iqr_bacteria_groups.(drug){j});   
                % join with Q1 table
                joined_bacteria_groups.(drug){j} = join(joined_bacteria_groups.(drug){j},q1_bacteria_groups.(drug){j});
                % join with Q3 table
                joined_bacteria_groups.(drug){j} = join(joined_bacteria_groups.(drug){j},q3_bacteria_groups.(drug){j});
            end
        end


        %~~~~~~~~~~~~~SYTO~~~~~~~~~~~~~%
        if syto_tf
            for i = drugs
                drug = i{1};
                for j = 1:size(iqr_syto_groups.(drug), 2)
                    % join average table and variance table 
                    joined_syto_groups.(drug){j} = join(center_syto_groups.(drug){j},iqr_syto_groups.(drug){j});   
                    % join with q1 table 
                    joined_syto_groups.(drug){j} = join(joined_syto_groups.(drug){j},q1_syto_groups.(drug){j});
                    % join with q3 table 
                    joined_syto_groups.(drug){j} = join(joined_syto_groups.(drug){j},q3_syto_groups.(drug){j});
                end
            end
        end 

        %~~~~~~~~~~~~~FM~~~~~~~~~~~~~%
        if FM_tf
            for i = drugs
                drug = i{1};
                for j = 1:size(iqr_FM_groups.(drug), 2)
                    % join average table and variance table 
                    joined_FM_groups.(drug){j} = join(center_FM_groups.(drug){j},iqr_FM_groups.(drug){j});    
                    % join with q1 table 
                    joined_FM_groups.(drug){j} = join(joined_FM_groups.(drug){j},q1_FM_groups.(drug){j});
                    % join with q3 table 
                    joined_FM_groups.(drug){j} = join(joined_FM_groups.(drug){j},q3_FM_groups.(drug){j});
                end
            end
        end 
    end 
end 


%% Add fluorescence together with bacteria and merge each drug into one table

all_data = cell2struct(cell(1,length(drugs)), drugs, 2);

has_all_img = true;

%merge files
for i = drugs
    drug = i{1};
    if cell_mode
        stop_point = length(grouped_no_nan_bacteria.(drug));
    else 
        stop_point = length(joined_bacteria_groups.(drug));
    end 

    % for following section, we used grouped_no_nan for cell workspace and
    % we use joined for img workspace. Boolean cell_mode is used to get from
    % one to the other. 
    for j = 1:stop_point    
        % get rows
        %find numeric cols for the avg files
        if cell_mode
            bact_row = grouped_no_nan_bacteria.(drug){j};
        else
            bact_row = joined_bacteria_groups.(drug){j};
        end
        numeric_avg_bacteria_cols = varfun(@isnumeric,bact_row,'OutputFormat', 'uniform');
        % if FM was selected earlier...
        if FM_tf
            % if cell workspace...
            if cell_mode
                FM_row = grouped_no_nan_FM_fluorescence.(drug){j};
                FM_row.f_NAME = FM_row.f_PARENT;
                FM_row.Properties.VariableNames{'f_NAME'} = 'NAME';
                % get duplicate values 
                % Unique values
                [~,idxu,idxc] = unique(FM_row.NAME);
                % count unique values
                [count, ~, idxcount] = histcounts(idxc,numel(idxu));
                % Where is greater than one occurence
                idxkeep = count(idxcount)>1;
                % Extract duplicate values
                dup_FM_row = FM_row(idxkeep,:);
                % now remove these from original FM row
                FM_row(idxkeep,:) = [];
                
                f_parents = unique(dup_FM_row.NAME)';
                for p = f_parents
                    f_name = p{1};
                    r_index = find(strcmp(dup_FM_row.NAME, f_name));
                    dup_inds = dup_FM_row(r_index,:);
                    % want to get a new row that has the mean values
                    numeric_cols = varfun(@isnumeric,dup_FM_row,'OutputFormat', 'uniform');
                    FM_data = dup_inds{:,numeric_cols};
                    non_numeric = dup_inds{1,~numeric_cols};
                    mean_FM_data = mean(FM_data,1);
                    averaged_row = [non_numeric mean_FM_data];
                    new_row = cell2table(cellstr(averaged_row),'VariableNames',dup_inds.Properties.VariableNames);
                    numeric_col_names = dup_inds.Properties.VariableNames(numeric_cols);
                    % convert from string to double for numeric values. 
                    for f = numeric_col_names
                        col_name = f{1};
                        new_row.(col_name) = str2double(new_row.(col_name));
                        
                    end
                    % add this new averaged row back
                    FM_row = [FM_row; new_row];
                end
            % else if image workspace 
            else
                FM_row = joined_FM_groups.(drug){j};
            end
            numeric_avg_FM_cols = varfun(@isnumeric,FM_row,'OutputFormat', 'uniform');
        end
        % if syto was selected earlier 
        if syto_tf
            % if cell workspace... 
            if cell_mode
                syto_row = grouped_no_nan_syto_fluorescence.(drug){j};
                syto_row.s_NAME = syto_row.s_PARENT;
                syto_row.Properties.VariableNames{'s_NAME'} = 'NAME';
                % get duplicate values 
                % Unique values
                [~,idxu,idxc] = unique(syto_row.NAME);
                % count unique values
                [count, ~, idxcount] = histcounts(idxc,numel(idxu));
                % Where is greater than one occurence
                idxkeep = count(idxcount)>1;
                % Extract duplicate values
                dup_syto_row = syto_row(idxkeep,:);
                % now remove these from original FM row
                syto_row(idxkeep,:) = [];
                
                s_parents = unique(dup_syto_row.NAME)';
                for p = s_parents
                    s_name = p{1};
                    r_index = find(strcmp(dup_syto_row.NAME, s_name));
                    dup_inds = dup_syto_row(r_index,:);
                    % want to get a new row that has the mean values
                    numeric_cols = varfun(@isnumeric,dup_syto_row,'OutputFormat', 'uniform');
                    syto_data = dup_inds{:,numeric_cols};
                    non_numeric = dup_inds{1,~numeric_cols};
                    mean_syto_data = mean(syto_data,1);
                    averaged_row = [non_numeric mean_syto_data];
                    new_row = cell2table(cellstr(averaged_row),'VariableNames',dup_inds.Properties.VariableNames);
                    numeric_col_names = dup_inds.Properties.VariableNames(numeric_cols);
                    % convert from string to double for numeric values. 
                    for f = numeric_col_names
                        col_name = f{1};
                        new_row.(col_name) = str2double(new_row.(col_name));
                        
                    end
                    % add this new averaged row back
                    syto_row = [syto_row; new_row];
                end 
                
            % else if image workspace... 
            else 
                syto_row = joined_syto_groups.(drug){j};
            end 
             numeric_avg_syto_cols = varfun(@isnumeric,syto_row,'OutputFormat', 'uniform');
        end
        % if cell workspace 
        if cell_mode
            id_col = bact_row.EXPERIMENT_fullname;
            img_col = bact_row.IMAGE_meta;
            drug_col = bact_row.EXPERIMENT;
        % else if image workspace 
        else
            id_col = bact_row.DRUG_id;
            img_col = bact_row.IMAGE_1;
            
            % get all image columns
            % old data does not have the all img so we need to create a
            % catch for that 
            try
                all_img_cols = contains(bact_row.Properties.VariableNames,"IMAGE");
            catch
                disp("no all img col")
                has_all_img = false; 
            end 
            all_q1_cols = contains(bact_row.Properties.VariableNames,"_q1");
            all_q3_cols = contains(bact_row.Properties.VariableNames,"_q3");
            all_iqr_cols = contains(bact_row.Properties.VariableNames,"_v");
            % don't want to get the repeats from q1 q3 and iqr
            if has_all_img 
                avg_img_cols = all_img_cols & ~all_q1_cols & ~all_q3_cols & ~all_iqr_cols;
                all_imgs = {bact_row{:,avg_img_cols}};
            end
            %%% could do something here for getting all images %%% 
            drug_col = drug;
        end

        
        drug_label = cell2table(cellstr(drug_col), 'VariableNames', {'DRUG'});
        id_label =  cell2table(cellstr(id_col), 'VariableNames', {'ID'});
        img_label = cell2table(cellstr(img_col), 'VariableNames', {'IMG'});
        if has_all_img
            all_img_label = table(all_imgs, 'VariableNames', {'ALL_IMG'});
        end
        %%%%%%%%%%% if cell workspace %%%%%%%%%%
        if cell_mode
            % if FM and syto were selected 
            if FM_tf & syto_tf
                sytoFM_row = join(bact_row,FM_row);
                all_row = join(sytoFM_row,syto_row);
                numeric_tot_cols = varfun(@isnumeric,all_row,'OutputFormat', 'uniform');
                kept_row = all_row(:,numeric_tot_cols);
                new_row = [drug_label id_label img_label kept_row];
            % if only FM was selected 
            elseif FM_tf
                all_row = join(bact_row,FM_row);
                numeric_tot_cols = varfun(@isnumeric,all_row,'OutputFormat', 'uniform');
                kept_row = all_row(:,numeric_tot_cols);
                new_row = [drug_label id_label img_label kept_row];
            % if only syto was selected
            elseif syto_tf
                all_row = join(bact_row,syto_row);
                numeric_tot_cols = varfun(@isnumeric,all_row,'OutputFormat', 'uniform');
                kept_row = all_row(:,numeric_tot_cols);
                new_row = [drug_label id_label img_label kept_row];
            % if neither syto nor FM were selected 
            else
                numeric_tot_cols = varfun(@isnumeric,bact_row,'OutputFormat', 'uniform');
                kept_row = bact_row(:,numeric_tot_cols);
                new_row = [drug_label id_label img_label kept_row];
            end
        %%%%%%%%% else if image %%%%%%%%%
         else 
            % if both syto and FM were selected
            if FM_tf & syto_tf
                new_row = [drug_label id_label img_label all_img_label bact_row(:,numeric_avg_bacteria_cols) FM_row(:,numeric_avg_FM_cols) syto_row(:,numeric_avg_syto_cols)];
            % if only FM was selected 
            elseif FM_tf
                new_row = [drug_label id_label img_label all_img_label bact_row(:,numeric_avg_bacteria_cols) FM_row(:,numeric_avg_FM_cols)];
            % if only syto was selected
            elseif syto_tf
                new_row = [drug_label id_label img_label all_img_label bact_row(:,numeric_avg_bacteria_cols) syto_row(:,numeric_avg_syto_cols)];
            % if neither syto nor FM were selected 
            else
                new_row = [drug_label id_label img_label all_img_label bact_row(:,numeric_avg_bacteria_cols)];
            end 
         end 
        all_data.(drug) = [all_data.(drug) ; new_row];
    end
end

clear num_bacteria_rows num_FM_rows num_syto_rows numeric_avg_bacteria_cols numeric_avg_FM_cols numeric_avg_syto_cols drug_label
clear new_row bact_row FM_row syto_row

%% Merge all tables into big final table

%create table
drug_tables = struct2cell(all_data);
final_data_table = vertcat(drug_tables{:});
if ~cell_mode
    final_data_table = removevars(final_data_table,config.removed_cols);
end 
%sort so averages and variances are next to each other
SortedNames = sort(final_data_table.Properties.VariableNames(4:end));
final_data_table = [final_data_table(:,1:3) final_data_table(:,SortedNames)];

%% Add dates column

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');

% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);

% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols); 

% Create date column 
% parsing of example image
% mm_03272019_drug52_set2_rep3_d__03_R3D_BaSiC
dates = extractAfter(final_data_table.IMG,"_"); % 03272019_drug52_set2_rep3_d__03_R3D_BaSiC
dates = extractBefore(dates,"_"); % 03272019
dates = strcat("D",dates); % D03272019
final_data_table.DATE = dates; 
final_data_table = [non_numeric final_data_table(:,end) numeric_data];

% get date indicies 
dates = cellstr(unique(dates))';
date_indexes = cell2struct(cell(1,length(dates)),dates,2);

for i = dates
    drug = i{1};
    date_indexes.(drug) = find(strcmp(final_data_table.DATE, drug)).';
end


%% Get batches for each 
%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');

% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);

% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols); 

% want to get batches from set, rep and pad


% parse through to get set_rep_pad
% parsing of example image
% mm_03272019_drug52_set2_rep3_d__03_R3D_BaSiC
batches = extractAfter(final_data_table.IMG,"_"); % 03272019_drug52_set2_rep3_d__03_R3D_BaSiC
batches = extractAfter(batches,"_"); % drug52_set2_rep3_d__03_R3D_BaSiC
batches = extractAfter(batches,"_"); % set2_rep3_d__03_R3D_BaSiC
sets = extractBefore(batches,"_"); % set2
batches = extractAfter(batches,"_"); % rep3_d__03_R3D_BaSiC
reps = extractBefore(batches,"_"); % rep3
batches = extractAfter(batches,"_"); % d__03_R3D_BaSiC
pads = extractBefore(batches,"_"); % d

batches = cellstr(strcat(sets, "_",reps, "_", pads)); %set2_rep3_d

set_rep = cellstr(strcat(sets, "_",reps)); %set2_rep3

% add batch column
final_data_table.BATCH = batches; %set2_rep3_d

% add rep column
final_data_table.REP = reps;

% add set_rep column
final_data_table.SET_REP = set_rep; %set2_rep3



% rearrange non-numeric and numeric information
final_data_table = [non_numeric final_data_table(:,end) ...
    final_data_table(:,end-1) final_data_table(:,end-2) numeric_data];

% get batch indicies 
batches = unique(batches)';
batch_indexes = cell2struct(cell(1,length(batches)),batches,2);

for i = batches
    drug = i{1};
    batch_indexes.(drug) = find(strcmp(final_data_table.BATCH, drug)).';
end

% get set_rep indicies
set_reps = unique(set_rep)';
set_rep_indexes = cell2struct(cell(1,length(set_reps)),set_reps,2);

for i = set_reps
    drug = i{1};
    set_rep_indexes.(drug) = find(strcmp(final_data_table.SET_REP, drug)).';
end

%% Add set column

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');

% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);

% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols); 

% parse through to get set_rep_pad
% parsing of example image
% mm_03272019_drug52_set2_rep3_d__03_R3D_BaSiC

% mm_07242019_x025_set1_rep3_a_02__5_R3D_BaSiC
no_mm = extractAfter(final_data_table.IMG,"_"); % 03272019_drug52_set2_rep3_d__03_R3D_BaSiC
% 07242019_x025_set1_rep3_a_02__5_R3D_BaSiC
no_date = extractAfter(no_mm,"_"); % drug52_set2_rep3_d__03_R3D_BaSiC
no_exp = extractAfter(no_date,"_"); % set2_rep3_d__03_R3D_BaSiC
sets = extractBefore(no_exp,"_"); % set2
% add set column
final_data_table.SET = sets; 
% rearrange numeric and non numeric information 
final_data_table = [non_numeric final_data_table(:,end) numeric_data];

% get set indicies 
sets = cellstr(unique(sets))';
set_indexes = cell2struct(cell(1,length(sets)),sets,2);

for i = sets
    drug = i{1};
    set_indexes.(drug) = find(strcmp(final_data_table.SET, drug)).';
end
    

%%  drug date column  
%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');

% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);

% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols); 

% add drug_date column
final_data_table.DRUG_DATE = cellstr(strcat(final_data_table.DRUG,"_",final_data_table.DATE));    
final_data_table = [non_numeric final_data_table(:,end) numeric_data];

% get drug_date indicies 
drug_dates = unique(final_data_table.DRUG_DATE)';
drug_date_indexes = cell2struct(cell(1,length(drug_dates)), drug_dates, 2);

for i = drug_dates
    drug = i{1};
    drug_date_indexes.(drug) = find(strcmp(final_data_table.DRUG_DATE, drug)).';
end


%% EXP column

% img name in style: mm_07242019_x025_set1_rep3_a_02__5_R3D_BaSiC
%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');

% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);

% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols); 

no_mm = extractAfter(final_data_table.IMG,"_"); % 07242019_x025_set1_rep3_a_02__5_R3D_BaSiC
no_date = extractAfter(no_mm,"_");% x025_set1_rep3_a_02__5_R3D_BaSiC

exp_col =  extractBefore(no_date,"_");% x025

% add exp column
final_data_table.EXP = exp_col;

% rearrange numeric and non numeric information 
final_data_table = [non_numeric final_data_table(:,end) numeric_data];

% get exp indicies
exps = cellstr(unique(exp_col))';
exp_indexes = cell2struct(cell(1,length(exps)),exps,2);

for i = exps
    drug = i{1};
    exp_indexes.(drug) = find(strcmp(final_data_table.EXP, drug)).';
end
    


%% normalization
%normalize final data table 
normalized_table = final_data_table;

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,normalized_table,'OutputFormat', 'uniform');

%find column names
numeric_final_col_names = normalized_table.Properties.VariableNames(numeric_final_data_cols);

%Final check and replace for NaN ***OR CAN USE OPTIONS IN PCA FOR DIFFERENT WAYS OF DEALING WITH NaN
normalized_table = replace_nan(normalized_table, numeric_final_data_cols);
% normalize data
normalized_table{:, numeric_final_data_cols} = normalize(normalized_table{:, numeric_final_data_cols},1,'range'); % normalize by dividing by largest value 

zero_mean_table = normalized_table; 
% 0 mean data
zero_mean_table{:, numeric_final_data_cols} = bsxfun(@minus,zero_mean_table{:, numeric_final_data_cols},mean(zero_mean_table{:, numeric_final_data_cols}));

% get numeric data from normalized table
ndata = zero_mean_table{:,numeric_final_data_cols}; 

drugs = unique(final_data_table.DRUG)';

%find row indexes for each drug
drug_indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

for i = drugs
    drug = i{1};
    drug_indexes.(drug) = find(strcmp(zero_mean_table.DRUG, drug)).';
end

% get unique images
imgs = unique(zero_mean_table.IMG)';

% get row indexes for each image
img_indexes = cell2struct(cell(1,length(imgs)), imgs, 2);

for i = imgs
    drug = i{1};
    img_indexes.(drug) = find(strcmp(zero_mean_table.IMG, drug)).';
end


batches = unique(final_data_table.BATCH)';
batches = cellstr(batches);
batch_indexes = cell2struct(cell(1,length(batches)),batches,2);

for i = batches
    drug = i{1};
    batch_indexes.(drug) = find(strcmp(final_data_table.BATCH, drug)).';
end

dates = unique(final_data_table.DATE)';
dates = cellstr(dates);
date_indexes = cell2struct(cell(1,length(dates)),dates,2);

for i = dates
    drug = i{1};
    date_indexes.(drug) = find(strcmp(final_data_table.DATE, drug)).';
end

% get id indicies
ids = unique(final_data_table.ID)';
ids = cellstr(ids);

% get row indicies for ids 
try 
id_indexes = cell2struct(cell(1,length(ids)),ids,2);
for i = ids
    drug = i{1};
    id_indexes.(drug) = find(strcmp(zero_mean_table.ID, drug)).';
end
catch
    disp("Unable to get id indexes--check for a . in id name")
end

% get a numerical list for each drug (eg if there are 3 Cer rows, will
% create a cell array including values Cer1 Cer2 Cer3). This is used in the
% creation of a clustergram. 
yvalues = (zero_mean_table.DRUG);
for i = drugs
    drug = i{1};
    for j = 1:length(drug_indexes.(drug))
        indx = drug_indexes.(drug)(j);
        yvalues(indx) = strcat(yvalues(indx),num2str(j));
    end
end

fprintf('Created final data table. Total time elapsed: %.2fs \n',toc)

 %% Make sure there are enough observations to plot PCA
if size(zero_mean_table,1) < 4 %must be 4 or more in order to have 3 Principal Components
    error(strcat('There are not enough points to calculate PCA! (', num2str(size(final_data_table,1)), ', must be 4 or greater). Make group size smaller! (Currently: ', num2str(config.group_size), ')'));
end

%% PCA Calculation
% Rows in the scores variable corresponds to each image/point inputed into the pca function
% Columns of the scores variable corresponds to the principal components (from first to last, decreasing order)
[coeff,scores,pcvars] = pca(ndata);
vari = cumsum(pcvars)./sum(pcvars);
coeff=coeff';
pcvari = pcvars/sum(pcvars);

numeric_final_col_names_space = numeric_final_col_names;
%Alters the column names (convert underscores (_) to spaces)
for i=1:size(numeric_final_col_names,2)
    numeric_final_col_names_space(1,i) = strrep(numeric_final_col_names(1,i),'_',' ');
end

% find all principal components and percentages for axes labels 
pc_labels = cell(1,length(pcvari));
for i = 1:length(pcvari)
    pc_labels(i) = cellstr(strcat("PC:",num2str(i)," ", (string(numeric_final_col_names_space(1,find(coeff(:,i) == max(coeff(:,i)))))), " ", (num2str(pcvari(i)*100)), "%"));
end

scores_table = horzcat(non_numeric,array2table(scores));

%% Create correlation plot 
vec = [100; 87.5; 75; 62.5; 50; 37.5; 25; 12.5; 0];
hex = ['#260000';'#610000';'#b00000';'#ff7575';'#feffff';'#87e2ff';'#00b0e9';'#006687';'#002b39'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 128;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
    
if config.output.show_correlation_plot
    % obtain correlation matrix
    corr_matrix = corrcoef(lda_data);

    % create heatmap
    figure();
    imagesc(corr_matrix);
    set(gca, 'XTick', 1:length(numeric_final_col_names_space)); % center x-axis ticks on bins
    set(gca, 'YTick', 1:length(numeric_final_col_names_space)); % center y-axis ticks on bins
    set(gca, 'XTickLabel', numeric_final_col_names_space); % set x-axis labels
    set(gca, 'YTickLabel', numeric_final_col_names_space); % set y-axis labels
    title('Correlation Matrix', 'FontSize', 14); % set title
    colorbar % create colorbar
    caxis([-1 1]) %set colorbar from -1 to 1
    ax = gca; % get current axes
    set(ax,'XTickLabelRotation',90) %rotate x axis labels 90 degrees

    % apply new colormap
    colormap(map)
end 

%% Clear or save

% clear these enormous variables 
clear transverse transverse_file transverse_no_nan

if config.output.save
    % get individual diretory name
    last_slash_pos = find(dname == '/', 1, 'last');
    folder_name = dname(last_slash_pos:end);
    if merge_tech_reps
        save_title = strcat(folder_name,"_merged_workspace.mat");
    else
        save_title = strcat(folder_name,"_workspace.mat");
    end 
    save(strcat(dname,save_title))
else 
    clear transverse_no_nan
end


fprintf('Finished. Total time elapsed: %.2fs \n',toc)