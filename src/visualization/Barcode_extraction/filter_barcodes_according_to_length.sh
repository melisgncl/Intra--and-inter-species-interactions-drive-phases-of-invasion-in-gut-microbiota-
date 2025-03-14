#!/bin/bash

# Check if at least 2 arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_fasta> <sequence_length1> [sequence_length2] [sequence_length3] ..."
    exit 1
fi

# Input arguments
INPUT_FASTA=$1
shift
SEQUENCE_LENGTHS=("$@")

# Output file to hold concatenated sequences with line numbers
OUTPUT_FILE="filtered_sequences_with_line_numbers.txt"

# Initialize counters for each sequence length
declare -A COUNT

# Process the FASTA file and filter reads with specified sequence lengths
awk -v seq_lengths="${SEQUENCE_LENGTHS[*]}" '
    BEGIN {
        split(seq_lengths, lengths_array);
        for (i in lengths_array) {
            lengths[lengths_array[i]] = 1;
            count[lengths_array[i]] = 0;
        }
        line_number = 0;
    }
    /^>/ {
        if (seq) {
            seq_len = length(seq);
            if (seq_len in lengths) {
                print seq "," line_number > "filtered_sequences_with_line_numbers.txt";
                count[seq_len]++;
            }
        }
        header = $0;
        seq = "";
        next;
    }
    {
        line_number++;
        seq = seq $0;
    }
    END {
        seq_len = length(seq);
        if (seq_len in lengths) {
            print seq "," line_number > "filtered_sequences_with_line_numbers.txt";
            count[seq_len]++;
        }
        for (len in lengths) {
            print "Number of reads with " len " characters: " count[len] > "/dev/stderr";
        }
    }
' "$INPUT_FASTA"

# Print summary
for SEQ_LEN in "${SEQUENCE_LENGTHS[@]}"; do
    echo "Number of reads with $SEQ_LEN characters: ${COUNT[$SEQ_LEN]}"
done

