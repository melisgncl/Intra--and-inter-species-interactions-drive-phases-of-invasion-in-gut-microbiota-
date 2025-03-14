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

ROOTDIR="/Users/melis/Dropbox/mouse_gut/github/"
mkdir -p "$ROOTDIR/reports/$Outdir/clustering/"

echo "Creating directories for clustering data..."
for cond in `cat condition2.list`; do
    for SAMPLE in `cat ${cond}.list`; do
    echo $cond $SAMPLE
    mkdir -p "$ROOTDIR/reports/$Outdir/clustering/$SAMPLE"
    done
done

echo "Plotting clusters and loess averages for a given treshold..."

##rm##
cp $ROOTDIR/reports/$Outdir/clustering_control/rm1/0.69/* $ROOTDIR/reports/$Outdir/clustering/rm1/
cp $ROOTDIR/reports/$Outdir/clustering_control/rm2/0.56/* $ROOTDIR/reports/$Outdir/clustering/rm2/
cp $ROOTDIR/reports/$Outdir/clustering_control/rm3/0.63/* $ROOTDIR/reports/$Outdir/clustering/rm3/
cp $ROOTDIR/reports/$Outdir/clustering_control/rm4/0.43/* $ROOTDIR/reports/$Outdir/clustering/rm4/

cp $ROOTDIR/reports/$Outdir/clustering_control/gf1/0.22/* $ROOTDIR/reports/$Outdir/clustering/gf1/
cp $ROOTDIR/reports/$Outdir/clustering_control/gf2/0.45/* $ROOTDIR/reports/$Outdir/clustering/gf2/
cp $ROOTDIR/reports/$Outdir/clustering_control/gf3/0.38/* $ROOTDIR/reports/$Outdir/clustering/gf3/
cp $ROOTDIR/reports/$Outdir/clustering_control/gf4/0.34/* $ROOTDIR/reports/$Outdir/clustering/gf4/



################
#### STEP 6 ####
################
echo "Generate dendrograms of clusters of chosen treshold  for each group..."
Rscript 6_cluster_of_clusters.R

################
#### STEP 7 ####
################
echo "Plot overlap index and simulate z scores within the cohort..."

Rscript 7_cluster_similarity.R
