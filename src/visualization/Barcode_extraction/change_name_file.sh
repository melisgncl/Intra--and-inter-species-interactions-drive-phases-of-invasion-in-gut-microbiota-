#!/bin/bash

# Check if directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Define the directory from the command-line argument
directory="$1"

# Loop through all relevant files in the specified directory
for file in "$directory"/*d*_S*barcode.txt; do
    # Extract the filename without the path
    filename=$(basename "$file")
    
    # Extract the numeric part following 'd'
    number=$(echo "$filename" | grep -oP 'd\K\d+')

    # Calculate the new number by adding 3
    new_number=$((number + 3))

    # Create the new file name
    new_file=$(echo "$filename" | sed -E "s/d${number}/t${new_number}/")

    # Rename the file (include the directory path)
    mv "$file" "$directory/$new_file"
    echo "Renamed: $filename -> $new_file"
done

