# MorphEUS 

Contains two main directories:

## segmentation_pipeline

Inside this directory is the code used to transform the .csv files outputted by MicrobeJ into a MATLAB workspace. This code filters out of focus and blurry cells using the transverse profile and calculates the median, Q1, Q3 and IQR for each replicate.


## classification_and_analysis

Inside this directory is the code used to run the classification trials as outlined in figure S3 as well as code to create the consensus KNN (cKNN) and heatmaps to visualize the results from the classification trials. 

In addition, code is included to perform the Kruskal–Wallis tests to determine signfiicance of variables based on using single-cell data and not replicate averages as seen in Figure 1. There is also code to create PCAs that were used in the earlier versions of the analysis pipeline (as seen in Figures 2 and S1). 

## Citation
If you find this work useful for your research, please cite our [paper](https://www.biorxiv.org/content/10.1101/2020.03.11.987545v1).