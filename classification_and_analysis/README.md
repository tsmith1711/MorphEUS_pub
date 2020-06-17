# MorphEUS Classification Trials and cKNNs

Performs a classification trial on a workspace created with MicrobeJ_segmentation

Currently set to use default workspaces used in the MophEUS paper, can be adjusted.

Run the script classification_trials.m to run a set of 70 classification trials. Results will be automatically saved into trial_results with the name the user selected along with a timestamp. A log containing all of the output produced during hte run will be saved into the logs folder. 

### Before you run

In order to run a classification trial, you need to download Hanchuan Peng's mRMR feature selection algorithm and compile it to work on your operating system. It can be found here:
[mRMR Feature Selection](https://www.mathworks.com/matlabcentral/fileexchange/14608-mrmr-feature-selection-using-mutual-information-computation?s_tid=prof_contriblnk). Place the resulting mRMR_0.9_compiled folder in the classification_trials directory.


### Workspaces for current presets

Joint profile: joint_workspace.mat

Timecourse applying: all_timecourse_doseresponse.mat

High dose: high_dose_workspace.mat

Low dose: low_dose_workspace.mat

Condition data: all_condition_dat_with_unt.mat

## Examining the data

To look through a set of classification trials, run cKNN_creator.m and select the desired workspace. 

To create boxplots using cell data created with MicrobeJ_segmentation, run boxplots_and_KW.m and selecte the desired workspace(s).

To create a PCA, run PCA_creator.m. 
