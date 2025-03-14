# Barcode Sequencing Data Processing
This folder contains scripts for processing barcode sequencing data.

## File Descriptions
*_condition.list: Lists the conditions (i.e., the mice to analyze). A separate list is created for each cohort.
*_sample.list: Lists the sequencing samples for analysis. Each cohort has its own list.
## Processing Steps
###Extract Barcodes
Run the script 1_extract_barcodes_new_with_blast2.sh to extract barcodes from the raw sequencing data.

###Error Correction and Time Point Separation

Run 2_run_Delcorrection.sh to correct errors within sequencing samples across all time points.
After error correction, separate the data by time points to ensure consistent barcode identities across the dataset.
Following these steps ensures accurate barcode identification and consistency across all time points.