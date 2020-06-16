%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MorphEUS classification trials %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script to perform classification trials and MorphEUS classifications
%%% Description: 
% prompts user if doing joint or individual profile, asks for other
% settings, asks user to select drugs to include in calculations, asks user
% to select drug(s) to apply, then performs MorphEUS classification trials
% saves result at end of maxrun number runs
% loads joint profile, takes 80 untreated from each 025x and 3x workspace,
% runs backwards select, picks out the best result (decided by highest
% percent success, if there is a tie, then with the lowest variable number)
% records results from what was pointing to what, then repeat the process
% from the top with a new random set of untreated
% will save result in trial_results and a log of all outputs in logs


%% clear workspace and set up paths 
clear

add_all_paths

rng('default')

% turn off TeX interpreter 
set(groot, 'defaultAxesTickLabelInterpreter', 'none')

%%  get user input for name of variable to save at end 
prompt = {'Enter a name for the variable information to be stored'};
dlgtitle = 'Save bayes_struct as:';
definput = {'bayes_struct'};
user_save_title = inputdlg(prompt,dlgtitle,[1 60],definput);

%% create a log
% get time string for a unique log name/variable name at the end 
time_string = strrep(datestr(clock),":","-");

logname = strcat("./logs/bayes log ",user_save_title, time_string,".txt");

diary(logname)
% for clarity of reading diary 
warning off 

%% some other settings 
with_randomized_labels = false;
% we are doing backwards select, so do feature reduction is true
do_feature_reduction = true;
% Turn off generation of images, off to not create (will be invisible with off)
create_cm = 'off'; 
show_pca = false; 

% settings for knn_helper
tit_add = "";
knn_type = 'median';
knn_with_apply = false;

% set this to know to save information after each trial 
bayesian = true;

% for new52 joint profile
use_default = false;

% names of workspaces to use
high_dose_wksp = "high_dose_workspace";
low_dose_wksp =  "low_dose_workspace";
joint_wksp = 'joint_workspace.mat';

%% joint or individual profile question
joint_profile_question = questdlg('Joint or individual profiles', ...
	'joint profiles?', ...
	'joint','individual','joint');
switch joint_profile_question
    case 'joint'
        do_joint_profile = true;
    case 'individual'
        do_joint_profile = false; 
end 

% if not joint profile, 3x or 025x
if ~do_joint_profile
    dose_question = questdlg('3x or 025x', ...
	'dose?', ...
	'3x','025x','3x');
    switch dose_question
        case '3x'
            x3_dose = true;
        case '025x'
            x3_dose = false; 
    end 
    
else % if do joint profiles, do flattened or separate?    
    flatten_question = questdlg('Flattened or separate?', ...
        '???',...
        'flattened','separate','flattened');
    switch flatten_question
        case 'flattened'
            flattened_joint = true;
        case 'separate'
            flattened_joint = false; 
    end
end 

% the usual 
variables_to_create = {'disc',  'include_DMSO', 'large_group','nn2',  'no_TVN','apply_condition_data','apply_timecourse_data','replicates_random','random_at_end','debug_mode'};
default_values =      [true,      false,          true,        false, false,    false,               false,          false,              false,     false];
disp("VARIABLE EXPLANATION")
disp("--")
disp("include_DMSO = include DMSO in TVN transformation")
disp("large group = categorize by larger groups instead of finer groups")
disp("nn2 = use second nearest neighbors as well")
disp("no_TVN = do not TVN transform")
disp("apply_condition_data = apply pH or col or but data (only works for individual 3x profiles)")
disp("apply_timecourse_data = apply 2h or 6h or 24h data (only works for individual profiles)")
disp("replicates_random = run this with random labels for each replicate")
disp("random_at_end = randomize the labels at the end, after doing backwards select")
disp("debug_mode = run in debug mode")
disp("--")

target_values = checkboxList("Select Settings", variables_to_create, default_values);

initialize_variables


% we've been using essentially this set so let's have it get
% suggested
suggestion = {'Mer','Amp','Ctax','INH','EMB','ETA','IMI','Van','Cyc','Del',...
'Lev','Mox','Clz','MIT','Olf','Kan','Amk','Cam','Cla','Dox','Gent',...
'Strep','Tet','Lin','Pre','CCCP','Cer','Mon','Nig','Thi','RifT','BDQ',... 
'RIF','THL','water','Untreated'};
%% do first run separate to not mess with parfor

% if debugging, only do 5 runs 
% num rows = number of untreated rows to keep 
% if debugging, do a small amount so it can run faster 
if debug_mode
    maxrun = 2; 
    num_rows = 32;
else
    maxrun = 70;
    num_rows = 80;
end 
runct = 1;

% don't want to just copy paste code again, this is the same content as the
% loop just over separately 
classification_loop_content 

%% run the rest of them 
% ah oh well we tried parfor is not working bby
for runct = 2:maxrun

    % since the exact same as before, just copied into a separate document 
    classification_loop_content 
   
end


%% auto save because this will take a while

if do_joint_profile
   disp("JOINT PROFILE")
   save_title = strcat("./trial_results/",user_save_title," joint ", time_string);
else
    disp("INDIVIDUAL PROFILE")
    if x3_dose
        disp("3x dose")
         save_title = strcat("./trial_results/",user_save_title," 3x ", time_string);
    else
        disp("025x dose")
         save_title = strcat("./trial_results/",user_save_title," 025x ", time_string);
    end
end

% want to indicate in the title if we had random replciates 
if replicates_random
   save_title = strcat(save_title, " replicates randomized"); 
end

if random_at_end
    save_title = strcat(save_title, " random at end"); 
    
end 

save(save_title,'bayes_struct','all_variable_sets')

disp(strcat("Saved bayes_struct as ", save_title))

diary off 