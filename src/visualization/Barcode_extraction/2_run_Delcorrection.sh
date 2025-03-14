#Adrian Serohijos 

###### POOL BARCODES ######
#	Since we have multiple time points per condition (mouse), we pool them to get a consistent list of clusters
#	For each condition (mouse) in list:
#	1. Append extracted barcodes from all time points to file (pooled-timepoints.txt)
#	2. cluster the pooled barcodes
###########################
for folder in `cat $1` ;
do
    mkdir -p $folder
    for COND in `cat ${folder}_condition.list`;
    do
        echo "Processing..." ${COND}
	mkdir -p $folder/$COND
	mkdir $folder/$COND/extracted_barcodes
	mkdir -p $folder/$COND/pooled_barcodes

        ##get extracted read from regex_bartender folder
        cp ../new_barcode_barcodeCounter2//$folder/extracted_barcodes/$COND*/*barcode.txt $folder/$COND/extracted_barcodes
        ###change d to t
	bash change_name_file.sh $folder/$COND/extracted_barcodes
        ###count the similar barcodes
	 ### Count the similar barcodes
        for barcode_file in $folder/$COND/extracted_barcodes/*barcode.txt; do
            output_file="${barcode_file%.txt}_count.txt"
            echo "Counting barcodes in $barcode_file"
            awk -F, '
            {
                # Remove the second column
                barcode = $1
                barcodes[barcode]++
            }
            END {
                # Print the barcodes with their counts
                for (barcode in barcodes) {
                    print barcode "\t" barcodes[barcode]
                }
            }
            ' "$barcode_file" > "$output_file"
        done

	rm $folder/$COND/extracted_barcodes/*barcode.txt

        echo "Running Error_Correction"
        python3 DeletionCorrection_edited3.py $folder/$COND/extracted_barcodes $folder/$COND/pooled_barcodes


    done
done
                 
