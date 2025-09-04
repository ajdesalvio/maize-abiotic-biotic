## Biotic and abiotic stress monitoring in maize using remote sensing and functional data analysis

This repository contains supplemental files and code walkthroughs to reproduce the analyses in our book chapter that explores the use of functional data analysis in maize biotic and abiotic stress monitoring with remote sensing.
Our book chapter can be found read at [this link](https://insert_link_here).

#### Case Study 1
- See the `Case Study 1` folder for the Python and R scripts used in the maize biotic stress (southern rust) case study.
The `Case Study 1 Data` folder contains files needed to run the scripts.
- To set up the python environment that was used in Case Study 1, please use the Case_Study_1.yml file to set up an Anaconda virtual environment.
- To set up the environment, first download Anaconda Navigator, and once installed, launch Anaconda Prompt. Navigate to the directory where the .yml file is saved using `cd path/to/file` and perform the installation as follows:
`conda env create -f environment.yml`
- Five python scripts are associated with Case Study 1 and their functions can be summarized as:
1) Create a shapefile from a CSV file. To do this, you'll need a CSV file that contains unique plot identifiers for each plot in the experiment. The coordinates to an "AB" line provides a way for the script to know which direction to orient the plots.
2) The shapefile created in the first step is next used to crop out individual plots within a large orthomosaic file (.tif).
3) GeoTIFFs produce cropped plots with large portions of black background, usually because the spatial orientation of the plot is retained after the cropping step. To remove the black background, this script identifies the region of interest that contains the plot and trims the black background. Currently this script only works for RGB TIFs.
4) Soil pixels can be segmented from plant pixels using the HUE or ExG indices, or an index provided by the user.
5) RGB vegetation indices can then be extracted from the segmented images. For each index, the mean and median value of the index is returned.

#### Case Study 2
See the `Case Study 2` folder for the Python and R scripts used in the maize abiotic stress (optimal vs. delayed planting) case study that uses Genomes to Fields data from Texas in 2021.
The `Case Study 2 Data` folder contains files needed to run the scripts.
The python environment setup for Case Study 2 can be performed similar to the procedure outlined in Case Study 1, except using the Case_Study_2.yml file.
