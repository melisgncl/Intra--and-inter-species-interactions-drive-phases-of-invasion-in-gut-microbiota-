# Import necessary libraries
import pandas as pd
import numpy as np
import time
from collections import Counter
import regex
import scipy.stats as sci_stats
import os
import sys

# Ensure the Levenshtein package is installed

import Levenshtein

# Define the get_deletion_neighborhood function
def get_deletion_neighborhood(stringer):
    return set([stringer] + [stringer[:x] + stringer[x+1:] for x in range(len(stringer))])

# Define the correct_bc_errors function
def correct_bc_errors(arr, min_counts_for_centroid=2, max_edits=3, poisson_error_rate=0.1):
    bc_counts = sorted(arr, key=lambda x: -1*x[1])
    error_correction_results = dict()
    deletion_dict = dict()
    seen_deletions = set()
    corrector = dict()

    for bc, count in bc_counts:
        del_net = get_deletion_neighborhood(bc)
        bc_parent = None
        dist_from_centroid = 0
        bc_parents = set()
        for d in del_net:
            if d in deletion_dict:
                bc_parents.add(deletion_dict[d])
        if bc_parents:
            if len(bc_parents) == 1:
                top_hit = bc_parents.pop()
            else:
                top_hit = sorted(bc_parents, key=lambda b: error_correction_results[b][2])[-1]
            if sci_stats.poisson.sf(count-1, error_correction_results[top_hit][2]*poisson_error_rate) > 0.05:
                bc_parent = top_hit
                dist_from_centroid = error_correction_results[top_hit][3]+1
                if dist_from_centroid <= max_edits:
                    while error_correction_results[top_hit][1] != 'centroid':
                        top_hit = error_correction_results[top_hit][1]
                    if Levenshtein.distance(top_hit, bc) <= max_edits:
                        corrector[bc] = top_hit
                    else:
                        corrector[bc] = 'excluded error'
                else:
                    corrector[bc] = 'excluded error'
        if not bc_parent:
            if count >= min_counts_for_centroid:
                bc_parent = 'centroid'
                corrector[bc] = bc
            else:
                bc_parent = 'excluded low count'
                corrector[bc] = 'excluded low count'
        if 'exclude' not in bc_parent:
            error_correction_results[bc] = [bc, bc_parent, count, dist_from_centroid]
            for d in del_net:
                if d not in seen_deletions:
                    deletion_dict[d] = bc
            seen_deletions.update(del_net)
    
    trues = [i[2] for i in error_correction_results.values() if i[1]=='centroid']
    errors = [i[2] for i in error_correction_results.values() if i[1]!='centroid']
    print(f'Dataset: {len(bc_counts)} bcs, {np.sum([i[1] for i in bc_counts])} reads.')
    print(f'"True" barcodes: {len(trues)}, totaling {sum(trues)} reads')
    print(f'Error barcodes: {len(errors)}, totaling {sum(errors)} reads')
    return corrector

# Define the error_correct function
def error_correct(df):
    td = df[df['BC'].apply(lambda bc: len(bc) > 5)].reset_index(drop=True)
    c = correct_bc_errors(np.array(td[['BC', 'Count']]))

    td['true_bc'] = td['BC'].apply(lambda b: c[b])
    td['exclusion_status'] = td['true_bc'].apply(lambda b: 'excluded' if b in ['excluded low count', 'excluded error'] else 'included')
    ttd = td[td['exclusion_status'] == 'included']
    bc_counts = ttd[['true_bc', 'Count']].groupby('true_bc').sum().reset_index().rename(columns={'true_bc': 'Sequence'})
    bc_counts['ClusterID'] = bc_counts.index + 1
    td['BCIDs'] = td.apply(lambda row: str(row.name + 1), axis=1)
    bcids = td[['true_bc', 'BCIDs']].groupby('true_bc')['BCIDs'].agg(lambda x: ','.join(x)).reset_index().rename(columns={'true_bc': 'Sequence'})
    return bc_counts.merge(bcids, on='Sequence', how='inner')[['ClusterID', 'Sequence', 'Count', 'BCIDs']], td

# Define the main function
def main(input_folder, output_folder):
    # Create a directory for results if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get all the input files from the input folder
    input_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.txt')]
    
    # Combine data from all input files and sum their counts
    combined_data = pd.DataFrame(columns=['BC', 'Count', 'Source'])
    for input_file in input_files:
        print(f"Processing file: {input_file}")
        data = pd.read_csv(input_file, sep='\t', header=None, names=['BC', 'Count'])
        data['Source'] = os.path.basename(input_file)  # Track the source file
        combined_data = pd.concat([combined_data, data], ignore_index=True)
    
    # Sum counts of the same barcodes
    combined_data = combined_data.groupby('BC').sum().reset_index()

    # Perform error correction on combined data
    print("Performing error correction on combined data...")
    startTime = time.time()
    corrected_counts, full_data = error_correct(combined_data)
    print('Execution time: ', time.time() - startTime, 'seconds')

    # Save corrected counts to a file
    corrected_counts.to_csv(os.path.join(output_folder, 'Combined_DeletionCluster_clusters.txt'), sep='\t', index=False)
    full_data.to_csv(os.path.join(output_folder, 'Combined_full_results.txt'), sep='\t', index=False)

    # Separate corrected barcodes back into their respective original files based on their counts
    for input_file in input_files:
        file_name = os.path.basename(input_file)
        original_data = pd.read_csv(input_file, sep='\t', header=None, names=['BC', 'Count'])
        merged_data = full_data.merge(original_data, left_on='BC', right_on='BC', how='left', suffixes=('_corrected', '_original'))
        output_file = os.path.join(output_folder, file_name + '_corrected_counts.txt')
        merged_data.to_csv(output_file, sep='\t', index=False)

    # Create a summary file with true barcode, true count, cluster ID, and time
    summary_data = pd.DataFrame(columns=['True_Barcode', 'True_Count', 'Cluster_ID', 'Time'])
    source_prefix = None
    for input_file in input_files:
        file_name = os.path.basename(input_file)
        # Extract time value from the filename (e.g., 't11' from 'GFM1t11_S12_R1_001.fastq.gz_passed_counts.txt')
        try:
            print(f"Extracting time value from file name: {file_name}")
            time_value = int(''.join(filter(str.isdigit, file_name.split('_')[0].split('t')[1])))
            print(f"Extracted time value: {time_value}")
        except IndexError:
            time_value = -1  # Set to -1 if no time value is found
            print(f"Failed to extract time value from file name: {file_name}, setting time value to -1")
        if source_prefix is None:
            source_prefix = file_name.split('t')[0]
        original_data = pd.read_csv(input_file, sep='\t', header=None, names=['BC', 'Count'])
        merged_data = full_data.merge(original_data, left_on='BC', right_on='BC', how='left', suffixes=('_corrected', '_original'))
        merged_data['Time'] = time_value
        grouped_data = merged_data.groupby(['true_bc', 'Time']).agg({'Count_original': 'sum'}).reset_index()
        grouped_data['Cluster_ID'] = grouped_data.groupby('true_bc').ngroup() + 1
        grouped_data.rename(columns={'true_bc': 'True_Barcode', 'Count_original': 'True_Count'}, inplace=True)
        summary_data = pd.concat([summary_data, grouped_data], ignore_index=True)

    # Save the summary data with time values
    summary_data.to_csv(os.path.join(output_folder, source_prefix + '.txt'), sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    main(input_folder, output_folder)

