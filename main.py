import pandas as pd
import glob, os, subprocess
import report, csv, time, gzip

from Bio import SeqIO
from collections import defaultdict


data_joe_path = '/Users/flavia/Documents/Concordia/C_KUZM1124/2nd_round/fastq/Sample_Chr4p_TSG'

# library original, replaced by targets_joe_unique after removing duplicates
# targets_joe = '/Users/flavia/Documents/Concordia/C_KUZM1124/2nd_round/fastq/Sample_Chr4p_TSG/targets_joe.fasta' #have duplicates

# library used to identify the unique sequences
targets_joe_unique = '/Users/flavia/Documents/Concordia/C_KUZM1124/2nd_round/fastq/Sample_Chr4p_TSG/targets_joe_unique.fasta'
# targets_joe_control = '/Users/flavia/Documents/Concordia/C_KUZM1124/2nd_round/fastq/Sample_Chr4p_TSG/targets_joe_only_control.fasta'


data_brittany_path = '/Users/flavia/Documents/Concordia/C_KUZM1124/2nd_round/fastq/Sample_MitoCarrier'
targets_brittany = '/Users/flavia/Documents/Concordia/C_KUZM1124/2nd_round/fastq/Sample_MitoCarrier/targets_brittany.fasta'

targets_DRs = '/Users/flavia/Documents/Concordia/C_KUZM1124/2nd_round/fastq/Sample_MitoCarrier/targets_DRs.fasta'


def verify_library(template_file_path, path_sequencing_data):
    """
    Verify if the library have sequences duplicates
    :param template_file_path: Path to the FASTA file containing gRNA sequences.
    :param path_sequencing_data: Path to the sequencing data directory.
    :return: list of gRNAS with duplicate sequences and create a new fasta removing the duplicates.
    """
    print("\n--- Verifying Library for Duplicate Sequences ---")
    # --- Initialize time tracking ---
    start_time = time.time()

    # --- Initialize variables ---
    sequences = defaultdict(list)  # Key: sequence, Value: list of IDs
    current_id = None
    current_sequence_lines = []
    total_sequences_count = 0

    # --- Parse the FASTA file ---
    with open(template_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # --- Check if the line is a header (starts with '>') ---
            if line.startswith('>'):
                if current_id and current_sequence_lines:
                    sequence = "".join(current_sequence_lines)
                    sequences[sequence].append(current_id)
                    total_sequences_count += 1

                current_id = line[1:].split()[0]  # Get the first part of the ID
                current_sequence_lines = []
            else:
                current_sequence_lines.append(line)

        # --- Process the last sequence in the file ---
        if current_id and current_sequence_lines:
            sequence = "".join(current_sequence_lines)
            sequences[sequence].append(current_id) # Add the last sequence to the dictionary
            total_sequences_count += 1

    # --- Prepare data for CSV and gather statistics ---
    duplicate_data = []
    sequences_with_duplicates_count = 0
    unique_sequences_count = 0

    for seq, ids in sequences.items():
        if len(ids) > 1:  # This sequence has duplicates (more than one ID)
            duplicate_data.append({
                'ids': ', '.join(ids),
                'id_count': len(ids),  # New: Count of unique IDs for this sequence
                'sequence': seq
            })
            sequences_with_duplicates_count += 1
        else:
            unique_sequences_count += 1

    # --- Write duplicate sequences to CSV ---
    with open(os.path.join(path_sequencing_data, 'library_duplicates.csv'), 'w', newline='') as out_f:
        # --- Define the column headers including the new 'id_count' ---
        fieldnames = ['ids', 'id_count', 'sequence']
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)

        writer.writeheader()
        if duplicate_data:
            writer.writerows(duplicate_data)
        else:
            # If no duplicates, the file will just contain the header
            pass

    # --- Create a fasta file with unique sequences ---
    unique_fasta_file = os.path.join(path_sequencing_data, 'unique_sequences.fasta')
    with open(unique_fasta_file, 'w') as fasta_out:
        for seq, ids in sequences.items():
            if len(ids) == 1:
                fasta_out.write(f">{ids[0]}\n{seq}\n")

    # --- Print Summary Statistics ---
    print("\n--- Sequence Analysis Summary ---")
    print(f"Total number of sequences in the file: {total_sequences_count}")
    print(f"Total number of unique sequences found: {unique_sequences_count}")
    print(
        f"Total number of sequences with duplicates (same sequence, different IDs): {sequences_with_duplicates_count}")
    if duplicate_data:
        print(f"Duplicate sequences (if any) written to: library_duplicates.csv")
    else:
        print("No duplicate sequences found with different IDs. The output CSV contains only headers.")

    elapsed_time = time.time() - start_time
    print(f"Finished: Library verification completed in {elapsed_time:.2f} seconds.\n")


def verify_qc(path_sequencing_data, fastqc_results_folder, multiqc_folder):
    """
    Verify the quality control of the sequences.
    """
    # Path to sequencing data
    seq_files = glob.glob(os.path.join(path_sequencing_data, '*.fastq.gz'))
    if not seq_files:
        seq_files = glob.glob(os.path.join(path_sequencing_data, '*.fastq'))

    # Define output path to the QC report files
    fastqc_results = os.path.join(path_sequencing_data, fastqc_results_folder)
    multiqc_report = os.path.join(path_sequencing_data, multiqc_folder)

    if not os.path.exists(fastqc_results):
        os.makedirs(fastqc_results)
    if not os.path.exists(multiqc_report):
        os.makedirs(multiqc_report)

    # Run FastQC on all fastq files
    for fastq_file in seq_files:
        print(f"Running FastQC on {fastq_file}")
        subprocess.run(['fastqc', fastq_file, '-o', fastqc_results], check=True)

    print(f"Running MultiQC on {fastqc_results}")
    subprocess.run(['multiqc', fastqc_results, '-o', multiqc_report], check=True)
    print("Finished QC Analysis")


def length_distribution(path_fasta_files):
    """
    Calculate the length distribution of sequences in FASTA files using seqkit.
    """
    fasta_files = glob.glob(os.path.join(path_fasta_files, '*.fasta'))
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in the specified directory.")

    for fasta_file in fasta_files:
        print(f"Processing {fasta_file}")
        # Calculate length distribution using seqkit
        subprocess.run(['seqkit', 'watch', '--fields', 'ReadLen', fasta_file, '-O', 'len_distribution.png'], check=True)

        try:
            # Run the first command: seqkit seq
            seqkit_seq = subprocess.run(
                ['seqkit', 'seq', '--min-len', '192', '--max-len', '192', fasta_file],
                check=True,
                capture_output=True,
                text=True
            )

            # Run the second command: seqkit stats
            seqkit_stats = subprocess.run(
                ['seqkit', 'stats'],
                input=seqkit_seq.stdout,
                check=True,
                text=True
            )

            print(seqkit_stats.stdout)  # Print the stats output
        except subprocess.CalledProcessError as e:
            print(f"Error processing {fasta_file}: {e}")


def remove_adapters(path_sequencing_data):
    """
    Remove adapters from FASTQ files using Cutadapt.
    """
    r1_files = glob.glob(os.path.join(path_sequencing_data, '*_R1_*.fastq.gz'))
    r2_files = glob.glob(os.path.join(path_sequencing_data, '*_R2_*.fastq.gz'))

    # Create a directory for cleaned files
    cleaned_dir = os.path.join(path_sequencing_data, 'cleaned_reads')
    if not os.path.exists(cleaned_dir):
        os.makedirs(cleaned_dir)

    if len(r1_files) != len(r2_files):
        raise ValueError("The number of R1 and R2 files do not match.")

    else:
        # Get a pair of R1 and R2 files
        for r1_file, r2_file in zip(r1_files, r2_files):
            base_name_r1 = os.path.splitext(os.path.splitext(os.path.basename(r1_file))[0])[0]
            base_name_r2 = os.path.splitext(os.path.splitext(os.path.basename(r2_file))[0])[0]
            output_file_r1 = os.path.join(cleaned_dir, f"{base_name_r1}.cleaned.fastq.gz")
            output_file_r2 = os.path.join(cleaned_dir, f"{base_name_r2}.cleaned.fastq.gz")
            print(f"Removing adapters from {base_name_r1} and {base_name_r2} to {output_file_r1} and {output_file_r2}")
            subprocess.run(['cutadapt', '-a', 'AATGATACGGCGACCACCGAGATCTACAC', '-A', 'CAAGCAGAAGACGGCATACGAGAT', '-o', output_file_r1, '-p', output_file_r2, r1_file, r2_file, '--cores=4'], check=True)


def merge_r1_r2(path_sequencing_data, output_folder, merge_tool):
    """
    Merge R1 and R2 reads into a single file.
    """
    r1_files = glob.glob(os.path.join(path_sequencing_data, 'raw_data', '*_R1_*.fastq.gz'))
    r2_files = glob.glob(os.path.join(path_sequencing_data, 'raw_data', '*_R2_*.fastq.gz'))

    # Create a directory for merged files
    merged_dir = os.path.join(path_sequencing_data, output_folder)
    if not os.path.exists(merged_dir):
        os.makedirs(merged_dir)

    if len(r1_files) != len(r2_files):
        raise ValueError("The number of R1 and R2 files do not match.")

    else:
        for r1_file in r1_files:
            r2_file = r1_file.replace("R1", "R2")
            if r2_file in r2_files:
                base_name = os.path.splitext(os.path.splitext(os.path.basename(r1_file))[0])[0].replace("_R1", "")
                output_prefix = os.path.join(merged_dir, base_name)
                if merge_tool == 'pear':
                    print(f"Merging {r1_file} and {r2_file}")
                    # --- Paramters for PEAR ---
                    # -q: quality score threshold for trimming the low quality part of a read
                    # -n: minimum length of the merged reads
                    # -t: minimum length of the reads after trimming
                    # -v: minimum overlap size (default is 10 bp)
                    # -f: forward read file
                    # -r: reverse read file
                    # -o: output prefix for the merged reads

                    # -- merged_v1 parameters
                    # subprocess.run(['pear', '-q', '25', '-t', '20', '-f', r1_file, '-r', r2_file, '-o', output_prefix], check=True)

                    # -- merged v2 parameters
                    # subprocess.run(['pear', '-q', '25', '-n', '140', '-t', '74', '-v', '5', '-f', r1_file, '-r', r2_file, '-o', output_prefix], check=True)

                    # -- merged v3 parameters
                    # subprocess.run(['pear', '-q', '20', '-v', '5', '-t', '20', '-f', r1_file, '-r', r2_file, '-o', output_prefix], check=True)

                    # -- merged v4 parameters
                    # subprocess.run(['pear', '-q', '20', '-v', '5', '-t', '96', '-f', r1_file, '-r', r2_file, '-o', output_prefix], check=True)

                    # -- merged v5 parameters
                    subprocess.run(['pear', '-q', '15', '-v', '10', '-m', '192', '-n', '190', '-p', '0.001', '-f', r1_file,
                                    '-r', r2_file, '-o', output_prefix], check=True)

                else:
                    # -- merge files with pandaseq
                    subprocess.run(['pandaseq', '-f', r1_file, '-r', r2_file, '-o', '5', '-O', '10', '-L', '192', '-l', '192', '-w', output_prefix+'.fasta'], check=True)


def clean_fastq_file(input_fastq, output_fastq):
    """
    Clean a FASTQ file by removing sequences with invalid characters (e.g., whitespace).
    """
    invalid_count = 0
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        try:
            for record in SeqIO.parse(infile, 'fastq'):
                if ' ' not in str(record.seq):  # Check for whitespace in the sequence
                    SeqIO.write(record, outfile, 'fastq')
                else:
                    print(f"Skipping invalid sequence in {input_fastq}: {record.id} (contains whitespace)")
                    invalid_count += 1
        except ValueError as e:
            if "Whitespace is not allowed in the sequence." in str(e):
                print(
                    f"Error: The FASTQ file {input_fastq} contains sequences with whitespace, which is invalid according to the FASTQ format. The file will be processed line by line to identify and skip these entries.")
                infile.seek(0)  # Go back to the beginning of the file
                record_lines = []
                for line in infile:
                    record_lines.append(line.strip())
                    if len(record_lines) == 4:
                        sequence = record_lines[1]
                        header = record_lines[0]
                        if ' ' not in sequence:
                            outfile.write('\n'.join(record_lines) + '\n')
                        else:
                            print(
                                f"Skipping invalid sequence in {input_fastq}: {header.split()[0]} (contains whitespace)")
                            invalid_count += 1
                        record_lines = []
            else:
                raise  # Re-raise other ValueErrors

    if invalid_count > 0:
        print(f"Skipped {invalid_count} invalid sequences from {input_fastq}.")


def convert_fastq_to_fasta(path_sequencing_data, merged_reads_folder):
    """
    Convert merged fastq files to fasta format and skip the sequences with invalid characters (e.g., whitespace).
    """
    invalid_count = 0
    merged_reads_path = os.path.join(path_sequencing_data, merged_reads_folder)
    fastq_files = glob.glob(os.path.join(merged_reads_path, '*.assembled.fastq'))
    fasta_dir = os.path.join(merged_reads_path, 'fasta')
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)

    for fastq_file in fastq_files:
        base_name = os.path.splitext(os.path.basename(fastq_file))[0]
        fasta_file = os.path.join(fasta_dir, f"{base_name}.fasta")
        print(f"Converting {fastq_file} to {fasta_file}")
        with open(fastq_file, 'r') as infile, open(fasta_file, 'w') as outfile:
            try:
                for record in SeqIO.parse(fastq_file, 'fastq'):
                    if ' ' not in str(record.seq):  # Check for whitespace in the sequence
                        SeqIO.write(record, outfile, 'fasta')
                    else:
                        print(f"Skipping invalid sequence in {fastq_file}: {record.id} (contains whitespace)")
                        invalid_count += 1
            except ValueError as e:
                if "Whitespace is not allowed in the sequence." in str(e):
                    print(
                        f"Error: The FASTQ file {fastq_file} contains sequences with whitespace, which is invalid according to the FASTQ format. The file will be processed line by line to identify and skip these entries.")
                    infile.seek(0)  # Go back to the beginning of the file
                    record_lines = []
                    for line in infile:
                        record_lines.append(line.strip())
                        if len(record_lines) == 4:
                            sequence = record_lines[1]
                            header = record_lines[0]
                            if ' ' not in sequence:
                                outfile.write('\n'.join(record_lines) + '\n')
                            else:
                                print(
                                    f"Skipping invalid sequence in {fastq_file}: {header.split()[0]} (contains whitespace)")
                                invalid_count += 1
                            record_lines = []
                else:
                    raise  # Re-raise other ValueErrors

        if invalid_count > 0:
            print(f"Skipped {invalid_count} invalid sequences from {fastq_file}.")


def filter_reads(path_fasta_files):
    """
    Filter reads with a lenght paramenter
    :param path_sequencing_data:
    :param merged_reads_folder:
    :return: fasta file with filtered reads
    """
    print("\n--- Filtering Reads by Length ---")
    count_filtered = 0
    total_sequences = 0
    fasta_files = glob.glob(os.path.join(path_fasta_files, '*assembled.fasta'))
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in the specified directory.")

    # Iterate through all FASTA files in the directory
    for fasta_file in fasta_files:
        print(f"Processing {fasta_file}")
        # Remove assembled from the file name and replace by filtered
        filtered_fasta_filename = os.path.basename(fasta_file).replace('assembled', 'filtered')

        # Create the filtered fasta file path
        filtered_fasta_file = os.path.join(path_fasta_files, filtered_fasta_filename)

        print(f"Filtering sequences in {filtered_fasta_filename} and saving to {filtered_fasta_file}")

        subprocess.run(['seqkit', 'seq', '-m', '192', '-M', '192', fasta_file, '-o', filtered_fasta_file], check=True)
        # subprocess.run(['seqkit', 'stats', filtered_fasta_file], check=True)

        print(f"Filtering complete.")


def count_reads(path_file_data):
    """
    Count the number of reads and calculate the length average in a FASTA file.
    """
    fasta_files = glob.glob(os.path.join(path_file_data, '*.fasta'))

    for fasta_file in fasta_files:
        print(f"Counting reads in {fasta_file}")
        with open(fasta_file, 'r') as infile:
            result = subprocess.run(['seqkit', 'stats', '-T', fasta_file],
                           capture_output=True, text=True, check=True)
            # Parse the output
            print(result.stdout)  # Print the stats output for debugging
            lines = result.stdout.strip().split('\n')
            if len(lines) > 1:
                stats = lines[1].split('\t')  # Second line contains the stats
                num_seqs = int(stats[3])  # Number of sequences
                avg_len = float(stats[6])  # Average length
                return num_seqs, avg_len
            else:
                raise ValueError("Unexpected seqkit output format.")


def get_drs_order(sequence, drs, regions):
    """
    Returns the order in which the DRs appear in the given sequence.

    Args:
        sequence (str): The DNA sequence to analyze.
        drs (list): List of DR names (e.g., ['DR4', 'DR3', 'DR1', 'DRWT']).
        regions (list): List of corresponding regions for the DRs.

    Returns:
        list: The order of DRs as they appear in the sequence.
    """
    dr_order = []
    for dr, region in zip(drs, regions):
        if region in sequence:
            dr_order.append(dr)
    return dr_order


def count_drs_order(merged_reads_folder):
    """
    Reads a FASTA file from the merged_reads_folder, calculates the order of DRs for each sequence,
    and prints the count of each unique order.

    Args:
        merged_reads_folder (str): Path to the folder containing the FASTA file.
    """
    # Define DRs and their corresponding regions
    drs = ['DRWT', 'DR4', 'DR3', 'DR1']
    regions = [
        "TAATTTCTACTCTTGTAGAT", # DRWT
        "TAATTTCTACTATTGTAGAT", # DR4
        "AAATTTCTACTCTAGTAGAT", # DR3
        "TAATTTCTACTGTCGTAGAT", # DR1
    ]

    # Find the FASTA file in the folder
    fasta_files = [f for f in os.listdir(merged_reads_folder) if f.endswith('.fasta')]
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in the specified directory.")
    fasta_file = os.path.join(merged_reads_folder, fasta_files[0])

    # Initialize a dictionary to count DR orders
    dr_order_counts = defaultdict(int)

    # Parse the FASTA file and calculate DR orders
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        dr_order = get_drs_order(sequence, drs, regions)
        dr_order_counts[tuple(dr_order)] += 1

    # Print the results
    print("DR Order Counts:")
    for order, count in dr_order_counts.items():
        print(f"{' -> '.join(order)}: {count} times")

    ''' Return a fasta file with sequence that contains all DRs '''
    output_fasta_file = os.path.join(merged_reads_folder, 'drs_order.fasta')
    with open(output_fasta_file, 'w') as fasta_out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)
            dr_order = get_drs_order(sequence, drs, regions)
            if len(dr_order) == len(drs):
                ''' Write the sequence to the output FASTA file if it contains all DRs '''
                fasta_out.write(f">{record.id}\n{sequence}\n")
    return dr_order_counts


def count_drs(merged_reads_folder, output_fasta_file):
    print(merged_reads_folder, output_fasta_file)
    drs = ['DRWT', 'DR4', 'DR3', 'DR1']

    regions = [
        "TAATTTCTACTATTGTAGAT",
        "AAATTTCTACTCTAGTAGAT",
        "TAATTTCTACTGTCGTAGAT",
        "TAATTTCTACTCTTGTAGAT",
    ]

    read_counts = {region: 0 for region in regions}
    total_reads_with_any_region = 0
    total_reads_with_all_regions = 0
    processed_reads = 0
    fastq_files = []
    fasta_files = []
    fastqgz_files = []

    fastq_files = glob.glob(os.path.join(merged_reads_folder, '*.assembled.fastq'))
    if not fastq_files:
        fasta_files = glob.glob(os.path.join(merged_reads_folder, '*.fasta'))
        if not fasta_files:
            fastqgz_files = glob.glob(os.path.join(merged_reads_folder, '*.fastq.gz'))
            if not fastqgz_files:
                raise FileNotFoundError("No files found in the specified directory.")

    if output_fasta_file is not None:
        fasta_output = open(output_fasta_file, 'w')

    if fastq_files:
        for fastq_file in fastq_files:
            print(f"Processing {fastq_file}")
            with open(fastq_file, 'r') as f:
                read_id = ""  # Initialize read_id for each file
                for i, line in enumerate(f):
                    if i % 4 == 0:  # This is the sequence ID line in FASTQ
                        # Remove the leading '@' and store the full ID
                        read_id = line.strip().lstrip('@')
                    elif i % 4 == 1:  # Read sequence line in FASTQ
                        processed_reads += 1
                        read_sequence = line.strip()
                        found_in_this_read = False
                        all_regions_found = True
                        for region in regions:
                            if region in read_sequence:
                                read_counts[region] += 1
                                found_in_this_read = True
                            else:
                                all_regions_found = False
                        if found_in_this_read:
                            total_reads_with_any_region += 1
                        if all_regions_found:
                            total_reads_with_all_regions += 1
                            if output_fasta_file is not None:
                                # Save the sequence to a FASTA file with the original read ID
                                fasta_output.write(f">{read_id}\n{read_sequence}\n")
    if fasta_files:
        for fasta_file in fasta_files:
            print(f"Processing {fasta_file}")
            with open(fasta_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    processed_reads += 1
                    read_sequence = str(record.seq)
                    read_id = record.id
                    found_in_this_read = False
                    all_regions_found = True
                    for region in regions:
                        if region in read_sequence:
                            read_counts[region] += 1
                            found_in_this_read = True
                        else:
                            all_regions_found = False
                    if found_in_this_read:
                        total_reads_with_any_region += 1
                    if all_regions_found:
                        total_reads_with_all_regions += 1
                        if output_fasta_file is not None:
                            # Save the sequence to a FASTA file with the original read ID
                            fasta_output.write(f">{read_id}\n{read_sequence}\n")
    if fastqgz_files:
        for fastqgz_file in fastqgz_files:
            print(f"Processing {fastqgz_file}")
            with gzip.open(fastqgz_file, 'rt') as f:
                read_id = ""  # Initialize read_id for each file
                for i, line in enumerate(f):
                    if i % 4 == 0:  # This is the sequence ID line in FASTQ
                        # Remove the leading '@' and store the full ID
                        read_id = line.strip().lstrip('@')
                    elif i % 4 == 1:  # Read sequence line in FASTQ
                        processed_reads += 1
                        read_sequence = line.strip()
                        found_in_this_read = False
                        all_regions_found = True
                        for region in regions:
                            if region in read_sequence:
                                read_counts[region] += 1
                                found_in_this_read = True
                            else:
                                all_regions_found = False
                        if found_in_this_read:
                            total_reads_with_any_region += 1
                        if all_regions_found:
                            total_reads_with_all_regions += 1
                            if output_fasta_file is not None:
                                # Save the sequence to a FASTA file with the original read ID
                                fasta_output.write(f">{read_id}\n{read_sequence}\n")

            print("Counts per region:")
            for dr, region in zip(drs, regions):
                print(f"  {dr} - {region}: {read_counts[region]} reads")
            print(f"\nTotal reads containing ANY of the regions: {total_reads_with_any_region}")
            print(f"Total reads containing ALL regions: {total_reads_with_all_regions}")
            print(f"Total reads processed: {processed_reads}")

    return drs, regions, read_counts, total_reads_with_any_region, total_reads_with_all_regions, processed_reads


def extract_regions(merged_reads_folder, output_dir):
    """
    Extract specific regions from sequences in a FASTA file and save them to separate FASTA files.

    :param input_fasta: Path to the input FASTA file.
    :param output_dir: Directory to save the output FASTA files.
    """
    # Define regions and their corresponding output files
    regions = {
        "crrna1": (26, 45),
        "crrna2": (66, 85),
        "crrna3": (106, 125),
        "crrna4": (146, 165),
        "crrna-all": (26, 165),  # Combined region for all crRNAs
    }
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Open output files for writing
    output_files = {name: open(f"{output_dir}/{name}.fasta", "w") for name in regions}
    print(f"Output files will be saved in {output_dir}.")

    # Convert fasta.gz to fasta if needed
    # fasta_files = glob.glob(os.path.join(merged_reads_folder, '*.fasta.gz'))
    # if fasta_files:
    #     print(f"Found {len(fasta_files)} gzipped FASTA files. Converting to FASTA...")
    #     for gz_file in fasta_files:
    #         with gzip.open(gz_file, 'rt') as f_in:
    #             fasta_content = f_in.read()
    #             fasta_file = gz_file.replace('.gz', '')
    #             with open(fasta_file, 'w') as f_out:
    #                 f_out.write(fasta_content)
    #         print(f"Converted {gz_file} to {fasta_file}.")
    #     # Update fasta_files to point to the uncompressed files

    fasta_files = glob.glob(os.path.join(merged_reads_folder, '*.fasta'))
    print(f"Found {len(fasta_files)} FASTA files in {merged_reads_folder}.")
    # Parse the input FASTA file
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            for name, (start, end) in regions.items():
                # Extract the subsequence (1-based indexing, so subtract 1 from start)
                subseq = record.seq[start - 1:end]
                # Write the subsequence to the corresponding output file
                output_files[name].write(f">{record.id}_{name}\n{subseq}\n")

    # Close all output files
    for file in output_files.values():
        file.close()

    print(f"Regions extracted and saved to {output_dir}.")


def count_cr_rnas(crrna_folder):
    """
    Count how many times each sequence appears in crRNA FASTA files and save counts to individual CSV files.

    :param crrna_folder: Path to the folder containing crRNA FASTA files.
    """

    # Find all FASTA files in the folder
    fasta_files = glob.glob(os.path.join(crrna_folder, 'crrna*.fasta'))
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in the specified directory.")

        # Process each FASTA file
    for fasta_file in fasta_files:
        sequence_counts = defaultdict(int)
        print(f"Processing {fasta_file}")

        # Count sequences in the FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)
            sequence_counts[sequence] += 1

        # Create a CSV file for the current FASTA file
        csv_file = os.path.splitext(fasta_file)[0] + '_counts.csv'
        with open(os.path.join(crrna_folder, csv_file), 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Sequence", "Count"])
            # Write the sequence counts decrescending order
            sequence_counts = dict(sorted(sequence_counts.items(), key=lambda item: item[1], reverse=True))
            for sequence, count in sequence_counts.items():
                writer.writerow([sequence, count])

        print(f"Counts saved to {csv_file}")


def find_sequences_in_fasta(crrna_folder, targets_fasta):
    """
    Find if sequences from a CSV file are present in a FASTA file and save the results.

    :param crrna_folder: Path to the CSV file containing sequences and counts.
    :param targets_fasta: Path to the FASTA file containing target sequences.
    :return output_csv: Path to save the output CSV file with sequence names and counts.
    """
    output_csv = 'matched_sequences.csv'
    cr_rna_counts_csv = glob.glob(os.path.join(crrna_folder, 'crrna-all_counts.csv'))[0]

    # Read sequences and counts from the CSV file
    csv_sequences = {}
    with open(cr_rna_counts_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            csv_sequences[row['Sequence']] = int(row['Count'])

    # Read sequences from the FASTA file
    fasta_sequences = {}
    for record in SeqIO.parse(targets_fasta, "fasta"):
        fasta_sequences[str(record.seq)] = record.id

    # Find matches and save results
    with open(os.path.join(crrna_folder, output_csv), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sequence_Name', 'Count'])  # Write header
        for sequence, count in csv_sequences.items():
            if sequence in fasta_sequences:
                writer.writerow([fasta_sequences[sequence], count])

    print(f"Matched sequences saved to {os.path.join(crrna_folder, output_csv)}.")


def run_alignment(path_sequencing_data, template_file, merged_reads_folder, alignment_tool):
    """
    Run Alignment
    """
    fasta_files = glob.glob(os.path.join(path_sequencing_data, merged_reads_folder, 'fasta', '*.filtered.fasta'))
    if not fasta_files:
        fasta_files = glob.glob(os.path.join(path_sequencing_data, merged_reads_folder, '*.fasta'))
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in the specified directory.")

    # Use the first FASTA file (or handle multiple files as needed)
    fasta_file = fasta_files[0]

    if alignment_tool == 'mmseqs2':
        """
        Run alignment using MMseqs.
        """

        # --- Run MMseqs commands ---
        # Create a database from the FASTA file and the template file
        subprocess.run(['mmseqs', 'createdb', fasta_file, 'merged_fasta_file_db'], check=True)
        subprocess.run(['mmseqs', 'createdb', template_file, 'template_db'], check=True)

        # Search the database against the template file
        subprocess.run(
            ['mmseqs', 'search', 'template_db', 'merged_fasta_file_db', 'AlignResult_mmseqs', 'tmp', '--search-type', '3', '--threads', '8'], check=True)

        # Convert the search results to a tab-separated format
        subprocess.run(
            ['mmseqs', 'convertalis', 'template_db', 'merged_fasta_file_db', 'AlignResult_mmseqs', 'AlignResult_mmseqs.tsv',
             '--format-output', 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'],
            check=True)

        # Delete DB temp files
        os.remove('merged_fasta_file_db.m8')

    else:
        script_dir = '/Users/flavia/PycharmProjects/C_KUZM1124'
        merged_db_name = os.path.join(script_dir, 'merged_fasta_file_db')
        align_result_file = os.path.join(script_dir, 'AlignResult_blastn.tsv')

        # --- Run alignment using Blastn ---
        # Create a database from the FASTA file and the template file
        subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl', '-out', merged_db_name], check=True)

        # Run blastn to search the database against the template file
        subprocess.run(['blastn', '-query', template_file, '-db', merged_db_name, '-out', align_result_file,
                        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
                        '-num_threads', '8'], check=True)


def filter_alignment_results(filter_pident, tool):
    """
    Filter the alignment results to keep only the best hits.
    """
    if tool == 'mmseqs':
        alignment_file = 'AlignResult_mmseqs.tsv'
    elif tool == 'blastn':
        alignment_file = 'AlignResult_blastn.tsv'
    else:
        raise ValueError("Unsupported alignment tool. Use 'mmseqs' or 'blastn'.")
    # Load the dataframe
    df = pd.read_csv(alignment_file, sep='\t', header=None,
                     names=["qseqid", "sseqid", "pident", "length", "mismatch",
                            "gapopen", "qstart", "qend", "sstart", "send",
                            "evalue", "bitscore"])

    # Print the number of entries in the filtered DataFrame
    print(f"Number of entries: {len(df)}")

    # Filter using string comparisons
    df_filtered = df[(df["pident"] >= filter_pident) & (df["length"] == 140)]
    # df_filtered = df[(df["pident"] >= filter_pident) & (df["length"] <= 142) & (df["length"] >= 138)]

    # Print the number of entries in the filtered DataFrame
    print(f"Number of entries after filtering: {len(df_filtered)}")

    # Save the filtered DataFrame
    df_filtered.to_csv(f'AlignResult_{tool}_filter_{str(filter_pident)}.tsv', sep='\t', index=False)

    # Calculate qseqid counts
    qseqid_counts = df_filtered[
        'qseqid'].value_counts().to_dict()  # This line was added to define and calculate qseqid_counts

    # Save the counts to a csv file
    with open(f'qseqid_counts_{str(filter_pident)}.txt', 'w') as f:
        for qseqid, count in qseqid_counts.items():
            f.write(f"{qseqid}\t{count}\n")

    with open(f'qseqid_counts_{str(filter_pident)}.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sequence_Name', 'Count'])  # Write header
        for qseqid, count in qseqid_counts.items():
            writer.writerow(qseqid, count)

    # Print the counts
    # print(qseqid_counts)
    print(f"Filtering complete. Results saved to {alignment_file} and qseqid_counts.txt.")


def fix_count_reads(duplicates_file, count_qseqid_file):
    """
    Fix the count of reads in the qseqid_counts file by checking the file with duplicates.
    :param duplicates: csv file with columns one with list of qseqids,
                      the second with the count of duplicates and the third with the sequence
    :param count_qseqid: columns qseqid and count
    :return: fixed_count_qseqid
    """

    """
        Fix the count of reads in the qseqid_counts file by checking the file with duplicates.

        :param duplicates_path: Path to a CSV file with columns:
                                - a column containing lists of qseqids (e.g., 'qseqids_list')
                                - a column with the count of duplicates (e.g., 'duplicate_count')
                                - a column with the sequence (e.g., 'sequence')
        :param count_qseqid_path: Path to a CSV file with columns:
                                  - 'qseqid'
                                  - 'count'
        :return: A pandas DataFrame with the fixed counts for qseqids.
        """

    try:
        # Load the duplicates file
        df_duplicates = pd.read_csv(duplicates_file)

        # Load the qseqid counts file
        # Using header=None, sep='\t', and comment='#' as per your previous fix.
        df_count_qseqid = pd.read_csv(count_qseqid_file, header=None, sep='\t', comment='#')
        df_count_qseqid = df_count_qseqid.rename(columns={0: 'qseqid', 1: 'original_count'})

    except FileNotFoundError as e:
        print(f"Error: One of the input files was not found: {e}")
        return None
    except pd.errors.EmptyDataError as e:
        print(f"Error: One of the input files is empty: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading files: {e}")
        return None

    # --- Preprocessing df_duplicates to map each individual ID to its original group ---
    # Create a unique identifier for each original duplicate group (e.g., the first ID in its 'ids' list)
    df_duplicates['group_id_original'] = df_duplicates['ids'].apply(lambda x: x.split(',')[0].strip())

    # 1. Split the 'ids' column by comma into a list
    df_duplicates['ids_list'] = df_duplicates['ids'].str.split(',')

    # 2. Explode the DataFrame so each individual ID in 'ids_list' gets its own row.
    # This will associate each 'cleaned_id' with its 'group_id_original'.
    df_duplicates_exploded = df_duplicates.explode('ids_list').copy()  # Use .copy() to avoid SettingWithCopyWarning

    # 3. Clean up the individual IDs (remove leading/trailing whitespace)
    df_duplicates_exploded['cleaned_id'] = df_duplicates_exploded['ids_list'].str.strip()

    # Select only the relevant columns for merging: cleaned_id and its original group_id
    duplicate_group_map = df_duplicates_exploded[['cleaned_id', 'group_id_original']].copy()

    # --- Merge count data with duplicate group information ---
    # Merge df_count_qseqid with the duplicate_group_map
    # This will add 'group_id_original' to df_count_qseqid for any qseqid that is part of a duplicate group.
    # qseqids not in a group will have NaN in 'group_id_original'.
    merged_df = pd.merge(df_count_qseqid,
                         duplicate_group_map,
                         left_on='qseqid',
                         right_on='cleaned_id',
                         how='left')

    # Drop the redundant 'cleaned_id' column from the merge
    merged_df = merged_df.drop(columns=['cleaned_id'])

    # --- Calculate the effective divisor for each group ---
    # For rows that belong to a group (group_id_original is not NaN):
    # Count how many members of each group are actually present in the qseqid_counts file.
    # This is done by counting non-NaN 'qseqid's within each 'group_id_original' group.
    # Use transform('count') to get the count for each group and broadcast it back to the original DataFrame size.
    merged_df['actual_group_members_present'] = merged_df.groupby('group_id_original')['qseqid'].transform('count')

    # Fill NaN values for 'actual_group_members_present' with 1 for non-grouped items
    # If group_id_original is NaN, it means it's not a duplicate group member, so its divisor is 1.
    merged_df['actual_group_members_present'] = merged_df['actual_group_members_present'].fillna(1).astype(int)

    # --- Calculate the fixed count based on the new logic ---
    merged_df['fixed_count'] = merged_df.apply(
        lambda row: row['original_count'] / row['actual_group_members_present'],
        axis=1
    )

    # Convert fixed_count to integer (if desired, results might be floats after division)
    # Be aware of potential data loss if division results in non-whole numbers.
    # If fixed counts are expected to be integers, you might want to round before converting.
    # For now, let's keep it as float to preserve precision from division.
    # If you need integer, you can add `.round().astype(int)`
    # merged_df['fixed_count'] = merged_df['fixed_count'].round().astype(int)

    # Clean up the DataFrame to return only 'qseqid' and 'fixed_count'
    fixed_count_qseqid = merged_df[['qseqid', 'fixed_count']].copy()

    # Create a new txt file to save the fixed counts
    output_file = f'fixed_{count_qseqid_file}'
    with open(output_file, 'w') as f:
        # Write header
        f.write("qseqid\tfixed_count\n")
        for index, row in fixed_count_qseqid.iterrows():
            # Format float to a reasonable number of decimal places if needed
            f.write(f"{row['qseqid']}\t{row['fixed_count']:.0f}\n")  # .0f for no decimal places, adjust as needed

    return fixed_count_qseqid


def find_missing_targets(template_path, filter_pident, alignment_csv_file_path):
    """
    Identifies targets present in a FASTA file but missing from an alignment CSV file.

    Args:
        template_path (str): Path to the input FASTA file containing all target sequences.
        alignment_csv_file_path (str): Path to the alignment TSV file.
        output_csv_path (str): Path to save the CSV file with missing target names.
        output_fasta_path (str): Path to save the FASTA file with missing target sequences.
    """
    print(f"Starting find_missing_targets analysis for FASTA: {template_path} and TSV: {alignment_csv_file_path}")

    output_csv_path = f'missing_targets_{filter_pident}.csv'
    output_fasta_path = f'missing_targets_{filter_pident}.fasta'
    # --- 1. Get all target names from the FASTA file ---
    fasta_target_names = set()
    try:
        for record in SeqIO.parse(template_path, "fasta"):
            fasta_target_names.add(record.id)
        print(f"Found {len(fasta_target_names)} unique targets in FASTA file.")
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {template_path}")
        return None
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return None

    # --- 2. Get unique qseqid (aligned target names) from the TSV file ---
    aligned_target_names = set()
    try:
        # print("Reading aligned target names from TSV file...")
        # Read the TSV file, specifying column names as per your description
        df = pd.read_csv(alignment_csv_file_path, sep=',', header=0, low_memory=False)
        # Get unique values from the 'Sequence_Name' column
        aligned_target_names = set(df["Sequence_Name"].unique())

        print(f"Found {len(aligned_target_names)} unique aligned targets in TSV file.")
    except FileNotFoundError:
        print(f"Error: TSV file not found at {alignment_csv_file_path}")
        return
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return

    # --- 3. Find targets missing from the alignment results ---
    missing_targets = fasta_target_names - aligned_target_names
    print(f"Found {len(missing_targets)} targets missing from the alignment results.")

    if not missing_targets:
        print("No missing targets found. No output files will be created.")
        return

    # --- 4. Save the missing target names to a CSV file ---
    try:
        # print(f"Saving missing target names to {output_csv_path}...")
        missing_df = pd.DataFrame(list(missing_targets), columns=["Missing_Target_Name"])
        missing_df.to_csv(output_csv_path, index=False)
        # print(f"Missing target names saved to {output_csv_path}")
    except Exception as e:
        print(f"Error saving CSV file: {e}")

    # --- 5. Create a subset of the FASTA file with the missing targets ---
    try:
        # print(f"Creating subset FASTA file for missing targets: {output_fasta_path}...")
        with open(output_fasta_path, "w") as output_handle:
            for record in SeqIO.parse(template_path, "fasta"):
                if record.id in missing_targets:
                    SeqIO.write(record, output_handle, "fasta")
        # print(f"Subset FASTA file saved to {output_fasta_path}")
    except Exception as e:
        print(f"Error creating subset FASTA file: {e}")

    return len(missing_targets)


if __name__ == '__main__':
    filter_pident = 100  # Default value for filter_pident

    # --- new pipeline ---
    # Joe's data processing and analysis
    # --- merge_tool = 'pear'  # or 'pandaseq' if you prefer
    verify_library(targets_joe_unique, data_joe_path)
    verify_qc(os.path.join(data_joe_path, 'raw_data'), 'raw_fastqc_results', 'raw_multiqc_report')
    merge_r1_r2(data_joe_path, output_folder='merged_reads_pandas', merge_tool='pandaseq')
    verify_qc(os.path.join(data_joe_path, 'merged_reads_pandas'), 'merged_fastqc_results', 'merged_multiqc_report')
    count_drs(os.path.join(data_joe_path, 'raw_data'), output_fasta_file='None')
    count_drs_order(os.path.join(data_joe_path, 'merged_reads_pandas'))
    extract_regions(os.path.join(data_joe_path, 'merged_reads_pandas'), os.path.join(data_joe_path, 'merged_reads_pandas', 'fasta'))
    count_cr_rnas(os.path.join(data_joe_path, 'merged_reads_pandas', 'fasta'))
    find_sequences_in_fasta(os.path.join(data_joe_path, 'merged_reads_pandas', 'fasta'), targets_joe_unique)
    find_missing_targets(targets_joe_unique, filter_pident, os.path.join(data_joe_path, 'merged_reads_pandas', 'fasta', "matched_sequences.csv"))
    report.generate_report(data_joe_path, filter_pident, 'merged_reads_pandas')

    # --- new pipeline ---
    # Brittany's data processing and analysis
    # --- merge_tool = 'pear'  # or 'pandaseq' if you prefer
    verify_library(targets_brittany, data_brittany_path)
    verify_qc(os.path.join(data_brittany_path, 'raw_data'), 'raw_fastqc_results', 'raw_multiqc_report')
    merge_r1_r2(data_brittany_path, output_folder='merged_reads_pandas', merge_tool='pandaseq')
    verify_qc(os.path.join(data_brittany_path, 'merged_reads_pandas'), 'merged_fastqc_results', 'merged_multiqc_report')
    count_drs(os.path.join(data_brittany_path, 'raw_data'), output_fasta_file='None')
    count_drs_order(os.path.join(data_brittany_path, 'merged_reads_pandas'))
    extract_regions(os.path.join(data_brittany_path, 'merged_reads_pandas'), os.path.join(data_brittany_path, 'merged_reads_pandas', 'fasta'))
    count_cr_rnas(os.path.join(data_brittany_path, 'merged_reads_pandas', 'fasta'))
    find_sequences_in_fasta(os.path.join(data_brittany_path, 'merged_reads_pandas', 'fasta'), targets_brittany)
    find_missing_targets(targets_brittany, filter_pident, "matched_sequences.csv")
    report.generate_report(data_brittany_path, filter_pident, 'merged_reads_pandas')


