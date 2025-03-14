#Adrian Serohijos 
# edited DGL

# where $1 folder

# where $2 is the data directory

# where $3 is the data type (ex: no_drug)


###### EXTRACT BARCODES ######
#   For a given sequencing sample from the list:
#   1. unpack all fastq files
#   2. Quality filter it
#   3.Extract with blast
#.  4.get the barcodes +2/-2 length
##############################

for folder in `cat $1` ;
do
	 base_folder=$(echo $folder | sed 's/[0-9]*//g')
         export NEW_DIR=$base_folder
	 echo $NEW_DIR

	 cd $folder
	 echo "in" $folder
	 mkdir processed_data2 
	 ls  *fastq.gz > sample.list

	for SAMPLE in `cat sample.list`;
	do
		echo "Processing " $SAMPLE
    		gunzip -d $SAMPLE

    		#pool reads
    		echo "Pooling barcodes from replicates ... "
    		cat *.fastq > processed_data2/pooled.fastq

    		gzip *.fastq
    
   		cd processed_data2
                echo "filter quality  and turn to fasta"
		python3 ../../fastq_to_fasta_quality_trim.py pooled.fastq 30

                # extract the barcodes
   		echo "Extracting barcodes."
                python3 ../../barcodeExtracter.py -templateSeq ../../referenceseq -inputFasta filtered_pooled.fasta -outputDir ${SAMPLE}

		echo "pooled sequences +2 adn -2"
		
		bash ../../filter_barcodes_according_to_length.sh ${SAMPLE}filtered_pooled_barcode.fasta 13 14 15 16 17

		mv filtered_sequences_with_line_numbers.txt ${SAMPLE}_barcode.txt


		mkdir /home/melis/mousebarcoding-master/mousebarcoding-master/Melis_2024_new_barcode_extraction/$NEW_DIR/extracted_barcodes/
		mkdir  /home/melis/mousebarcoding-master/mousebarcoding-master/Melis_2024_new_barcode_extraction/$NEW_DIR/extracted_barcodes/$SAMPLE
                
                  
                mv *barcode.txt  /home/melis/mousebarcoding-master/mousebarcoding-master/Melis_2024_new_barcode_extraction/$NEW_DIR/extracted_barcodes/$SAMPLE/
                  
		rm * 

		cd ../
    
	done
	cd ../
done


