#!/usr/bin/env python

## Python script
## Louis Gauthier, August 2020

#######################    READ ME!   #####################################

##### This python script calculates the correlation matrix for a time series of filtered barcodes

#######################    ^^^^^^^^   #####################################

#args[1] input barcode trajectory file
#args[2] output path
#args[3] sample name (ex: SM1)

# ex usage: python src/visualization/5_clustering/2_apply_clustering.py data/clustering/GFM1-4/GFM2_filtered.csv data/clustering/GFM1-4/ GFM2

import csv
import sys
import math
import numpy as np
import numpy.ma as ma
import pandas as pd

################################################################################

## 1. read sample data
## INPUT: ROW VECTORS OF FREQUENCY TRAJECTORIES
sample_data = np.genfromtxt(sys.argv[1], delimiter=',',skip_header=1)

sample_path = sys.argv[2]
sample_name = sys.argv[3]

################################################################################

## 2. calculate distance matrix (pearson correlation)
## DISTANCE: 1 - (pearson correlation coefficient between two trajectories)
df = pd.DataFrame(sample_data[:,1:np.size(sample_data,1)-1])

print(len(df))

## 	  must exclude time points where either frequency = 0 to calculate the correlation coefficient
##	  thus we set all 0's to NaN
df[df == 0] = np.nan

log_data = np.log10(df)

##	  we mask NaN cells from calculation (handled by corrcoef)
## 	  WARNING: corrcoef is not guaranteed to return x in [-1,1] because of floating-point rounding 
corr_mat = ma.corrcoef(ma.masked_invalid(log_data))

##	  compute distance matrix
dist_mat = np.matrix(np.ones((len(df), len(df))) - corr_mat)

dist_mat[dist_mat < 0] = 0

##	  save distance matrix
pd.DataFrame(dist_mat).to_csv("%s/%s_dist.csv" % (sample_path, sample_name))




