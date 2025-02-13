ROOTDIR="/Users/melis/Dropbox/mouse_gut/github"

echo "Creating directories for clustering data..."
for taxa_cond in `cat taxa.list`;
do
    for taxa in `cat ${taxa_cond}.list`;
    do
        for clone_cond in `cat clone.list`;
        do
          for clone in `cat ${clone_cond}.list`;
           do
               for cutoff in `cat  cutoff.list`; 
               do
                   echo $taxa_cond $taxa $clone_cond $clone $cutoff ${taxa_cond}_${clone_cond}
                   mkdir -p $ROOTDIR/reports/barcodeCounter_Delcorrect/${taxa_cond}_${clone_cond}
                   mkdir -p $ROOTDIR/reports/barcodeCounter_Delcorrect/mixing_index/${taxa_cond}_${clone_cond}/${taxa}_${clone}
               done
           done
      done
  done
done







echo "Coclustering within and between cohorts and calculating mixing index...."
for taxa_cond in `cat taxa.list`;
do
    for taxa in `cat ${taxa_cond}.list`;
    do
        for clone_cond in `cat clone.list`;
        do
            for clone in `cat ${clone_cond}.list`;
            do
                for cutoff in `cat  cutoff.list`; 
                do
                    
                    sample="${taxa}_${clone}"
                    echo $taxa_cond $taxa $clone_cond $clone $cutoff $sample 
                    Rscript 2_coclustering_statistic_with_all_cutoff.R $ROOTDIR/data/16S/taxa_long_format/${taxa}_16S_taxa.tsv  $ROOTDIR/reports/barcodeCounter_Delcorrect/clustering_control/$clone/$cutoff/${clone}_clustered_loess_log10.csv $sample ${taxa_cond}_${clone_cond} $cutoff
                done
            done
       done
   done
done





echo "Ploting mixing index... "
Rscript 3_plot_mixing_index.R