import glob, os, subprocess

import pandas as pd
from Bio import SeqIO


data_joe_path = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_Chr4p_TSG'
targets_joe = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_Chr4p_TSG/targets_joe.fasta'

data_brittany_path = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_MitoCarrier'
targets_brittany = '/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_MitoCarrier/targets_brittany.fasta'


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


def count_reads(path_sequencing_data, merged_reads_folder):
    """
    Count the number of reads and calculate the length average in a FASTA file.
    """
    total_length = 0
    num_sequences = 0
    fasta_dir = os.path.join(path_sequencing_data, merged_reads_folder, 'fasta')
    fasta_file = glob.glob(os.path.join(fasta_dir, '*.fasta'))[0]
    # fasta_file = glob.glob(os.path.join('/Users/flavia/Documents/Concordia/C_KUZM1024/2nd_round/fastq/Sample_MitoCarrier', '*.fasta'))[0]

    with open(fasta_file, 'r') as infile:
        for record in SeqIO.parse(infile, 'fasta'):
            total_length += len(record.seq)
            num_sequences += 1
    if num_sequences > 0:
        average_length = total_length / num_sequences
        print(f'Average length of sequences in fasta: {average_length:.2f}')
    print(f'Number of sequences in fasta: {num_sequences}')


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

if __name__ == '__main__':
    # Joe's data processing and analysis
    verify_qc(data_joe_path, 'raw_fastqc_results', 'raw_multiqc_report')
    # remove_adapters(data_brittany_path) # Data does not contain adapters
    merge_r1_r2(data_joe_path)
    verify_qc(os.path.join(data_joe_path, 'merged_reads'), 'merged_fastqc_results', 'merged_multiqc_report')
    convert_fastq_to_fasta(data_joe_path, 'merged_reads')
    run_alignment(data_joe_path, targets_joe, 'merged_reads')
    filter_alignment_results()

    # Brittany's data processing and analysis
    verify_qc(os.path.join(data_brittany_path, 'raw_data'), 'raw_fastqc_results', 'raw_multiqc_report')
    # remove_adapters(data_brittany_path) # Data does not contain adapters
    merge_r1_r2(data_brittany_path)
    verify_qc(os.path.join(data_brittany_path, 'merged_reads'), 'merged_fastqc_results', 'merged_multiqc_report')
    convert_fastq_to_fasta(data_brittany_path, 'merged_reads')
    count_reads(data_brittany_path, 'merged_reads')
    run_alignment(data_brittany_path, targets_brittany, 'merged_reads')
    filter_alignment_results()