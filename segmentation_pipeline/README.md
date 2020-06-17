# MorphEUS segmentation 

Image segmentation process with data from MicrobeJ.

*See the [full documentation](https://tufts.box.com/s/070ctjfuuk8uuwo2sbngb3bw09mkbn13) for running this program in Box. (Somewhat outdated as of 6/20/19)*

### How to run

This script is run after using MocrbeJ to segment a 3 channel image (one brightfield, two flourescent) and exporting their output as csv files. The process of using MicrobeJ to segment an image is more closely detailed in the documentation linked above.

Once you have your results from segmentation, create a folder with a name that is meaningful to the data contained.  Inside this folder, create a folder called 'data', where you will place your csvs named with the following naming scheme:

all_bacteria.csv 
all_feature1.csv
all_feature2.csv
all_transverse.csv

When running the script MicrobeJ_segmentation, you will be prompted to select the parent folder. During this time, select your folder that contains the 'data' folder. Do NOT select the data folder itself. The script will automatically save a workspace with the name of the folder followed by '_merged_workspace.mat' inside the folder. 

Optional: Put a config file in the folder (the most up to date one is located in this repository, you can just copy and paste that one). If no config file is placed, by default the one in this repository will be used.

Optional: To combine several files for larger data sets, you can add additional sets with names such as all_bacteria_2.csv all_feature1_2.csv, etc. as long as they are all placed in the 'data' folder.

### The workspace

The most important thing in the workspace created will be saved in final_data_table. This table contians the results from filtering out out of focus cells and calculating the median, q1, q3 and iqr for each replicate. 

