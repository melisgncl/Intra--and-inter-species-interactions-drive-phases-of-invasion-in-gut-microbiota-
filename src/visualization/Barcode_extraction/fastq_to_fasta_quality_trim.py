import subprocess
import sys
import gzip
from Bio import SeqIO

def count_reads(file_path, file_format):
    """
    Count the number of reads in a given file (FASTQ or FASTA).
    """
    try:
        with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path, "r") as handle:
            return sum(1 for _ in SeqIO.parse(handle, file_format))
    except Exception as e:
        print(f"Error counting reads in {file_path}. Error: {e}")
        return 0

def run_fastp(input_fastq, quality_threshold, output_fastq):
    """
    Run fastp to filter reads based on quality threshold.
    """
    command = [
        "fastp",
        "-i", input_fastq,
        "-o", output_fastq,
        "-q", str(quality_threshold)
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Filtered {input_fastq} and saved to {output_fastq} with quality threshold {quality_threshold}")
    except subprocess.CalledProcessError:
        print("fastp failed. Please check the input file and try again.")
        sys.exit(1)

def convert_fastq_to_fasta(input_fastq, output_fasta):
    """
    Convert filtered FASTQ to FASTA.
    """
    try:
        with gzip.open(input_fastq, "rt") if input_fastq.endswith(".gz") else open(input_fastq, "r") as handle:
            sequences = SeqIO.parse(handle, "fastq")
            with open(output_fasta, "w") as fasta_handle:
                SeqIO.write(sequences, fasta_handle, "fasta")
        print(f"Converted {input_fastq} to {output_fasta}")
    except Exception as e:
        print(f"Failed to convert {input_fastq} to FASTA format. Error: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fastq> <quality_threshold>")
        sys.exit(1)

    input_fastq = sys.argv[1]
    quality_threshold = int(sys.argv[2])
    output_fastq = f"filtered_{input_fastq.split('/')[-1].replace('.fastq.gz', '').replace('.fastq', '')}.fastq.gz"
    output_fasta = f"filtered_{input_fastq.split('/')[-1].replace('.fastq.gz', '').replace('.fastq', '')}.fasta"

    # Count reads before filtering
    initial_read_count = count_reads(input_fastq, "fastq")

    run_fastp(input_fastq, quality_threshold, output_fastq)

    # Count reads after filtering
    filtered_read_count = count_reads(output_fastq, "fastq")

    convert_fastq_to_fasta(output_fastq, output_fasta)

    # Count reads in final FASTA
    final_fasta_read_count = count_reads(output_fasta, "fasta")

    # Print summary
    print("Summary:")
    print(f"Input FASTQ: {input_fastq}")
    print(f"Quality threshold: {quality_threshold}")
    print(f"Initial read count: {initial_read_count}")
    print(f"Filtered read count: {filtered_read_count}")
    print(f"Final FASTA read count: {final_fasta_read_count}")
    print(f"Filtered FASTQ: {output_fastq}")
    print(f"Output FASTA: {output_fasta}")

if __name__ == "__main__":
    main()

