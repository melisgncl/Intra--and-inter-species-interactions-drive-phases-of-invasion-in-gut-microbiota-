#!/usr/bin/env python

## Python script
## Louis Gauthier, August 2020

import csv
import sys
import math
import numpy as np
import matplotlib
matplotlib.use('agg')
import pylab
import scipy
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
import numpy.ma as ma
import pandas as pd
import collections

def flatten(iterable):
	for el in iterable:
		if isinstance(el, collections.Iterable) and not isinstance(el, str): 
			yield from flatten(el)
		else:
			yield el

def check_symmetric(a, rtol=1e-05, atol=1e-08):
	return np.allclose(a, a.T, rtol=rtol, atol=atol)

sample_dist = sys.argv[1]
sample_path = sys.argv[2]
sample_name = sys.argv[3]
threshold = float(sys.argv[4])

matplotlib.rcParams['lines.linewidth'] = 0.5

## READ DISTANCE MATRIX

dist_mat = np.genfromtxt(sys.argv[1], delimiter=',',skip_header=1)
dist_mat = dist_mat[:,1:]

# verify matrix symmetry with reasonable tolerance for floating-point error
#print(check_symmetric(dist_mat))

# convert vector-form distance matrix to square-form
# (squareform is dumb and will throw an error so we ignore checks)
square_dist_mat = ssd.squareform(dist_mat,checks=False)

################################################################################

# 3. calculate linkage matrix

# AVERAGE METHOD (UPGMA)
linkage_matrix = sch.linkage(square_dist_mat,method="average",optimal_ordering=True)

#print(linkage_matrix)

################################################################################

TICK_LABEL_SIZE=12
AXIS_LABEL_SIZE=12

# Initialize figure of clustering summary
pylab.rc("font", family="serif", size=TICK_LABEL_SIZE)

fig = pylab.figure(figsize=(4, 4))

# Dendrogram
ax = fig.add_axes([0,0.7,0.7,0.2], frame_on=True)

sch.set_link_color_palette(['hotpink','tan','teal','darkmagenta','springgreen','sienna','darkturquoise','darkkhaki','violet','darkorchid','crimson','darkorange','forestgreen','royalblue'])

# clustering threshold
#ct = 0.5*max(linkage_matrix[:,2])
ct = threshold

dendrogram = sch.dendrogram(linkage_matrix, count_sort='descending',ax=ax, color_threshold=ct, get_leaves=True, no_labels=True,above_threshold_color="gray")

ax.axhline(ct, linestyle="--", color="black")

ax.set_ylabel("Correlation", fontsize=0.7*AXIS_LABEL_SIZE)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
ax.set_yticklabels(["1.00", "0.75", "0.50", "0.25", "0.00", "-0.25", "-0.50"], fontsize=0.75*TICK_LABEL_SIZE)


ax.spines["right"].set_color('none')
ax.spines["top"].set_color('none')
ax.spines["bottom"].set_color('none')
ax.yaxis.set_ticks_position("left")

pylab.setp(ax.get_xticklines(), visible=False)
pylab.setp(ax.get_xticklabels(), visible=False)
ax.set_title(sample_name, y=1.2, fontweight="bold", fontsize=AXIS_LABEL_SIZE)

# Permute indices of distance matrix to match the dendrogram
for i in range(len(dist_mat)):
	dist_mat[:, i] = dist_mat[dendrogram["leaves"], i]
for i in range(len(dist_mat)):
	dist_mat[i, :] = dist_mat[i, dendrogram["leaves"]]

# Correlation matrix
ax = fig.add_axes([0,0,0.7,0.7])
im = ax.imshow(np.ones((len(dist_mat), len(dist_mat))) - dist_mat, cmap="bwr", vmin=-1, vmax=1)        
ax.set_xlabel("Lineage trajectories", fontsize=0.75*AXIS_LABEL_SIZE)
ax.set_ylabel("Lineage trajectories", fontsize=0.75*AXIS_LABEL_SIZE)
pylab.setp(ax.get_xticklabels(), visible=False)
pylab.setp(ax.get_yticklabels(), visible=False)
pylab.setp(ax.get_xticklines(), visible=False)
pylab.setp(ax.get_yticklines(), visible=False)

#ax.text(-0.2, 0.5, "side label", va="center", fontweight="bold", fontsize=0.8*AXIS_LABEL_SIZE, rotation=90, transform=ax.transAxes)

# Colorbar for correlation matrix
ax = fig.add_axes([0.7,0.1,0.2,0.5])
cbar = pylab.colorbar(im, cax=ax)
cbar.set_label("Correlation", fontsize=0.75*AXIS_LABEL_SIZE, labelpad=12, rotation=270)
ax.set_aspect(15)

fig.savefig("%s/%s_dendromap_average.pdf" % (sample_path,sample_name), bbox_inches="tight")
fig.savefig("%s/%s_dendromap_average.png" % (sample_path,sample_name), bbox_inches="tight", dpi=300)

clusters = sch.fcluster(linkage_matrix, t=ct, criterion='distance')

#print(clusters)
#print(flatten(clusters))
#print(dendrogram['leaves'])

rows = zip(flatten(clusters),[x+1 for x in dendrogram['leaves']])

################################################################################

with open("%s/clusters_%s_average" % (sample_path,sample_name), "w") as f:
	writer = csv.writer(f)
	for row in rows:
		writer.writerow(row)


## CLUSTERING DIAGNOSTICS ################################################################################

#matplotlib.pyplot.figure(figsize=(10,10))

#matplotlib.pyplot.scatter(corr_mat[:,0],corr_mat[:,1],c=clusters)

#matplotlib.pyplot.savefig("diagnostics.pdf", bbox_inches="tight")