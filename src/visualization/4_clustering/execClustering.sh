#!/bin/bash

###give outdir
# Function to display usage instructions
usage() {
  echo "Usage: $0 -o <output_directory>"
  exit 1
}

# Parse command-line arguments
while getopts ":o:" opt; do
  case ${opt} in
    o )
      Outdir=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

# Ensure Outdir is provided
if [ -z "${Outdir}" ]; then
  usage
fi

ROOTDIR="/Users/melis/Dropbox/mouse_gut/github"
PROC_SAMPLES="$ROOTDIR/data/$Outdir"
DATA_CLUST="$ROOTDIR/data/clustering_${Outdir}"

echo "Creating directories for clustering data..."
for cond in `cat condition2.list`; do
    for SAMPLE in `cat ${cond}.list`; do
        mkdir -p "$ROOTDIR/reports/$Outdir/clustering_control/"
        for cutoff in `cat cutoff.list`; do
            echo $cond $SAMPLE
            mkdir -p "$DATA_CLUST/$SAMPLE"
            mkdir -p "$ROOTDIR/reports/$Outdir/clustering_control/$SAMPLE"
            mkdir -p "$ROOTDIR/reports/$Outdir/clustering_control/$SAMPLE/$cutoff"
        done
    done
done

################
#### STEP 1 ####
################
echo "Filtering trajectories for clustering..."

# RM samples
Rscript 1_filter_data.R $PROC_SAMPLES/M1/Sample_M1_clustering.txt $DATA_CLUST/rm1/rm1_filtered.csv 0 12 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/M2/Sample_M2_clustering.txt $DATA_CLUST/rm2/rm2_filtered.csv 0 12 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/M3/Sample_M3_clustering.txt $DATA_CLUST/rm3/rm3_filtered.csv 0 12 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/M4/Sample_M4_clustering.txt $DATA_CLUST/rm4/rm4_filtered.csv 0 12 5e-5

# Germ-free samples
Rscript 1_filter_data.R $PROC_SAMPLES/GFM1/Sample_GFM1_clustering.txt $DATA_CLUST/gf1/gf1_filtered.csv 0 13 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/GFM2/Sample_GFM2_clustering.txt $DATA_CLUST/gf2/gf2_filtered.csv 0 13 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/GFM3/Sample_GFM3_clustering.txt $DATA_CLUST/gf3/gf3_filtered.csv 0 14 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/GFM4/Sample_GFM4_clustering.txt $DATA_CLUST/gf4/gf4_filtered.csv 0 13 5e-5

################
#### STEP 2 ####
################
echo "Computing correlation matrix..."
for cond in `cat condition2.list`; do
    for SAMPLE in `cat ${cond}.list`; do
        echo Sample $SAMPLE
        python 2_apply_clustering.py $DATA_CLUST/$SAMPLE/${SAMPLE}_filtered.csv $DATA_CLUST/$SAMPLE $SAMPLE
    done
done

################
#### STEP 3 ####
################
# See 3_dendrogram_upgma.py for details
# but remember that the last argument is the threshold for flattening the clusters
# this threshold is set manually by the user based on inspection of the dendrogram/heatmap
# to choose a threshold run a range between 0.1-0.9 with step size 0.02 and check smallest distance between loess then choose a cut off
# and plot clusters and loess averages
echo "Generating dendrogram and flattening clusters..."
for cond in `cat condition2.list`; do
    for SAMPLE in `cat ${cond}.list`; do
        for cutoff in `cat cutoff.list`; do 
            echo $cond $SAMPLE $cutoff
            python3 3_dendrogram_upgma.py $DATA_CLUST/$SAMPLE/${SAMPLE}_dist.csv $ROOTDIR/reports/$Outdir/clustering_control/$SAMPLE/$cutoff $SAMPLE $cutoff
            echo "Plotting clusters and loess averages..."
            Rscript 4_plot_clusters_loess_new_c.R $ROOTDIR/reports/$Outdir/clustering_control/$SAMPLE/$cutoff/clusters_${SAMPLE}_average $DATA_CLUST/$SAMPLE/${SAMPLE}_filtered.csv $ROOTDIR/reports/$Outdir/clustering_control/$SAMPLE/$cutoff/$SAMPLE $cond
        done
    done
done

################
#### STEP 4 ####
################
#### quantify the clustering to choose for best threshold to define clone dynamics
echo "Generate plot of hierarchical cluster quantification for each group..."
Rscript 5_h_clustering_quantification.R 
