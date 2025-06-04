import pandas as pd
import glob, os, subprocess
import report, csv

from Bio import SeqIO
from collections import defaultdict


data_joe_path = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_Chr4p_TSG'
targets_joe = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_Chr4p_TSG/targets_joe.fasta'

data_brittany_path = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_MitoCarrier'
targets_brittany = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_MitoCarrier/targets_brittany.fasta'


def verify_library(template_file):
    """
    Verify if the library have sequences duplicates
    :param template_file:
    :return: list of gRNAS with duplicate sequences
    """
    sequences = defaultdict(list)  # Key: sequence, Value: list of IDs
    current_id = None
    current_sequence_lines = []
    total_sequences_count = 0

    # --- Parse the FASTA file ---
    with open(template_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if current_id and current_sequence_lines:
                    sequence = "".join(current_sequence_lines)
                    sequences[sequence].append(current_id)
                    total_sequences_count += 1

                current_id = line[1:].split()[0]  # Get the first part of the ID
                current_sequence_lines = []
            else:
                current_sequence_lines.append(line)

        # Process the last sequence in the file
        if current_id and current_sequence_lines:
            sequence = "".join(current_sequence_lines)
            sequences[sequence].append(current_id)
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
    with open('library_duplicates.csv', 'w', newline='') as out_f:
        # Define the column headers including the new 'id_count'
        fieldnames = ['ids', 'id_count', 'sequence']
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)

        writer.writeheader()
        if duplicate_data:
            writer.writerows(duplicate_data)
        else:
            # If no duplicates, the file will just contain the header
            pass

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


def merge_r1_r2(path_sequencing_data):
    """
    Merge R1 and R2 reads into a single file.
    """
    r1_files = glob.glob(os.path.join(path_sequencing_data, '*_R1_*.fastq.gz'))
    r2_files = glob.glob(os.path.join(path_sequencing_data, '*_R2_*.fastq.gz'))

    # Create a directory for merged files
    merged_dir = os.path.join(path_sequencing_data, 'merged_reads')
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
                print(f"Merging {r1_file} and {r2_file}")
                # parameters for PEAR -q: quality threshold, -t: minimum length of the merged reads, -f: forward read file, -r: reverse read file, -o: output directory
                # the default minimun overlap size is set to 10 bp
                subprocess.run(['pear', '-q', '25', '-t', '20', '-f', r1_file, '-r', r2_file, '-o', output_prefix], check=True)


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


def count_reads(path_file_data):
    """
    Count the number of reads and calculate the length average in a FASTA file.
    """
    total_length = 0
    num_sequences = 0
    average_length = 0
    fasta_file = glob.glob(os.path.join(path_file_data, '*.fasta'))[0]

    with open(fasta_file, 'r') as infile:
        for record in SeqIO.parse(infile, 'fasta'):
            total_length += len(record.seq)
            num_sequences += 1
    if num_sequences > 0:
        average_length = total_length / num_sequences
        # print(f'Average length of sequences in fasta: {average_length:.2f}')
    # print(f'Number of sequences in fasta: {num_sequences}')
    return num_sequences, average_length


def run_alignment(path_sequencing_data, template_file, merged_reads_folder):
    """
    Run alignment using MMseqs.
    """
    fasta_files = glob.glob(os.path.join(path_sequencing_data, merged_reads_folder, 'fasta', '*.fasta'))
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in the specified directory.")

    # Use the first FASTA file (or handle multiple files as needed)
    fasta_file = fasta_files[0]

    # Run MMseqs commands
    subprocess.run(['mmseqs', 'createdb', fasta_file, 'merged_fasta_file_db'], check=True)
    subprocess.run(['mmseqs', 'createdb', template_file, 'template_db'], check=True)
    subprocess.run(
        ['mmseqs', 'search', 'template_db', 'merged_fasta_file_db', 'AlignResult_mmseqs', 'tmp', '--search-type', '3',
         '--threads', '4'], check=True)
    subprocess.run(
        ['mmseqs', 'convertalis', 'template_db', 'merged_fasta_file_db', 'AlignResult_mmseqs', 'AlignResult_mmseqs.tsv',
         '--format-output', 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'],
        check=True)


def filter_alignment_results():
    """
    Filter the alignment results to keep only the best hits.
    """
    # Set filter parameter
    filter_pident = 100  # Default value for filter_pident

    # Load the dataframe
    df = pd.read_csv('AlignResult_mmseqs.tsv', sep='\t', header=None,
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
    df_filtered.to_csv(f'AlignResult_mmseqs_filter_{str(filter_pident)}.tsv', sep='\t', index=False)

    # Calculate qseqid counts
    qseqid_counts = df_filtered[
        'qseqid'].value_counts().to_dict()  # This line was added to define and calculate qseqid_counts

    # Save the counts to a file
    with open(f'qseqid_counts_{str(filter_pident)}.txt', 'w') as f:
        for qseqid, count in qseqid_counts.items():
            f.write(f"{qseqid}\t{count}\n")

    # Print the counts
    # print(qseqid_counts)
    print("Filtering complete. Results saved to AlignResult_mmseqs_filter.tsv and qseqid_counts.txt.")


def find_missing_targets(template_path, alignment_tsv_file_path, output_csv_path, output_fasta_path):
    """
    Identifies targets present in a FASTA file but missing from an alignment TSV file.

    Args:
        template_path (str): Path to the input FASTA file containing all target sequences.
        alignment_tsv_file_path (str): Path to the alignment TSV file.
        output_csv_path (str): Path to save the CSV file with missing target names.
        output_fasta_path (str): Path to save the FASTA file with missing target sequences.
    """
    print(f"Starting find_missing_targets analysis for FASTA: {template_path} and TSV: {alignment_tsv_file_path}")

    # --- 1. Get all target names from the FASTA file ---
    fasta_target_names = set()
    try:
        for record in SeqIO.parse(template_path, "fasta"):
            fasta_target_names.add(record.id)
        # print(f"Found {len(fasta_target_names)} unique targets in FASTA file.")
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {template_path}")
        return
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return

    # --- 2. Get unique qseqid (aligned target names) from the TSV file ---
    aligned_target_names = set()
    try:
        # print("Reading aligned target names from TSV file...")
        # Read the TSV file, specifying column names as per your description
        df = pd.read_csv(alignment_tsv_file_path, sep='\t', header=None, low_memory=False,
                         names=["qseqid", "sseqid", "pident", "length", "mismatch",
                                "gapopen", "qstart", "qend", "sstart", "send",
                                "evalue", "bitscore"])
        # Get unique values from the 'qseqid' column
        aligned_target_names = set(df["qseqid"].unique())
        # print(f"Found {len(aligned_target_names)} unique aligned targets in TSV file.")
    except FileNotFoundError:
        print(f"Error: TSV file not found at {alignment_tsv_file_path}")
        return
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return

    # --- 3. Find targets missing from the alignment results ---
    missing_targets = fasta_target_names - aligned_target_names
    # print(f"Found {len(missing_targets)} targets missing from the alignment results.")

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
    # Joe's data processing and analysis
    verify_library(targets_joe)
    verify_qc(data_joe_path, 'raw_fastqc_results', 'raw_multiqc_report')
    ## remove_adapters(data_joe_path) # Data does not contain adapters
    merge_r1_r2(data_joe_path)
    verify_qc(os.path.join(data_joe_path, 'merged_reads'), 'merged_fastqc_results', 'merged_multiqc_report')
    convert_fastq_to_fasta(data_joe_path, 'merged_reads')
    length_distribution(os.path.join(data_joe_path, 'merged_reads', 'fasta'))
    run_alignment(data_joe_path, targets_joe, 'merged_reads')
    filter_alignment_results()
    find_missing_targets(targets_brittany, '', 'missing_targets.csv', 'missing_targets.fasta')
    report.generate_report(data_brittany_path, 'merged_reads', '')

    # Brittany's data processing and analysis
    verify_qc(os.path.join(data_brittany_path, 'raw_data'), 'raw_fastqc_results', 'raw_multiqc_report')
    ## remove_adapters(data_brittany_path) # Data does not contain adapters
    merge_r1_r2(data_brittany_path)
    verify_qc(os.path.join(data_brittany_path, 'merged_reads'), 'merged_fastqc_results', 'merged_multiqc_report')
    convert_fastq_to_fasta(data_brittany_path, 'merged_reads')
    length_distribution(os.path.join(data_brittany_path, 'merged_reads', 'fasta'))
    count_reads(os.path.join(data_brittany_path, 'merged_reads', 'fasta'))
    run_alignment(data_brittany_path, targets_brittany, 'merged_reads')
    filter_alignment_results()
    find_missing_targets(targets_brittany, 'AlignResult_mmseqs_filter_97.tsv', 'missing_targets.csv', 'missing_targets.fasta')
    report.generate_report(data_brittany_path, 'merged_reads', 'qseqid_counts_97.txt')