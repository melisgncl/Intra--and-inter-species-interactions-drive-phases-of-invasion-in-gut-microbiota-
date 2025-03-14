# Quantifying the intra- and inter-species community interaction in a microbiome by dynamic covariance mapping

This is the repository archiving data and scripts for reproducing results presented  "Quantifying the intra- and inter-species community interaction in a microbiome by dynamic covariance mapping"

## Developed by the <http://www.serohijoslab.org/>

![](https://github.com/melisgncl/Intra--and-inter-species-interactions-drive-phases-of-invasion-in-gut-microbiota-/blob/main/reports/Readme_Figures/DCM_overwiew.jpg?raw=true)

## Dynamic Covariance Mapping (DCM) Framework

The Dynamic Covariance Mapping (DCM) introduces a parameter-free methodology for estimating the community interaction matrix directly from abundance time-series data of microbial community members. This framework is designed to address the complexity of microbial interactions in their natural environments where multiple species interactions and intra-species diversity significantly influence the community dynamics.

######  

To estimate the community interaction matrix and analyze stability changes across different phases using Dynamic Covariance Mapping (DCM), simply execute the following script:src/general_DCM.R


## Setup and Configuration

### Dependencies and Libraries: 
Before beginning the analysis, ensure all necessary libraries are loaded. This can be done by  running the script located at:
 src/visualization/0_config/0_config.R
 

## Importing Data
	
### Barcode Clustering Data Import and Reshaping: 
To start analyzing your data, first import and reshape the barcode clustering data into a convenient format using scripts in:
  src/visualization/1_intersection/*


## Analysis and Visualization
	
### Barcode Dynamics
Visualize barcode dynamics using Muller-style area plots and log-line plots. Scripts for these visualizations are found at:
 src/visualization/2_dynamics/1_plotDynamics.R

### Barcode Diversity
Calculate the diversity of barcodes for all samples and plot the results with the following scripts:
 Calculate Diversity: src/visualization/3_diversity/1_calculateDiversity.R
 Plot Diversity: src/visualization/3_diversity/2_plotDiversity.R

### Clone Dynamics
Determine the dynamics of clones through analysis scripts located in:
 src/visualization/4_clustering/
	
### 16S rRNA Analysis
16S rRNA Gene Sequencing Data Analysis: Analyze and visualize 16S rRNA gene sequencing data to study bacterial compositions. Scripts for these analyses are available at:
 src/visualization/5_16S/

### Co-clustering Analysis
Co-clustering Community and Clone Dynamics: Examine the co-clustering of community dynamics with clone dynamics using scripts in:
 src/visualization/6_coclustering/

### Dynamical Covariance Mapping (DCM) Analysis
DCM Analysis: Perform Differential Condition Matrix analyses to further understand the conditions affecting the microbiome. Relevant scripts are located at:
 src/visualization/7_DCM/