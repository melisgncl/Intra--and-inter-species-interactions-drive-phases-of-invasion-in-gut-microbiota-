# Barcode Composition

This folder contains files for the composition of barcodes in samples.

## Directory Structure
### raw 
The raw folder contains barcodes extracted from FASTQ files with blast. The data is organized in tabular format as follows:

Barcode_sequence: Shows the extracted barcode sequence.\
Remaining Columns: Indicate the time points with their respective counts.

### processed 
The processed folder contains putative barcodes identified after processing with Deletion_Correct. The data is organized in tabular format as follows:

Center: Shows the putative barcode sequence.\
Remaining Columns: Indicate the time points with their respective counts.
