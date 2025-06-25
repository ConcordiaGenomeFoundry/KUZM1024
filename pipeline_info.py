from reportlab.platypus import Paragraph, Spacer, PageBreak
from reportlab.lib.units import inch


def generate_report(story, styles, filter_pident):
    story.append(Paragraph("<b>Computational Pipeline</b>", styles['ReportSectionHeading']))
    story.append(Paragraph(
        "The following steps outline the workflow to process raw illumina sequencing reads and verify the reads against a library database. ",
        styles['ReportBodyText']))
    story.append(Spacer(1, 0.3 * inch))

    pipeline_title = [
        "1. Raw Reads Quality Verification",
        "2. Reads Merging",
        "3. Merged Reads Quality Verification"
        "4. Verification of DRs in Merged Reads",
        "5. Extract Internal Sequence from Reads",
        "6. Process Internal Sequence",
        "7. Verify Library Database",
    ]

    pipeline_steps = [
        "Initial quality control of raw sequencing reads was performed using FastQC to generate individual read quality reports. These reports were then aggregated and summarized with MultiQC for a comprehensive overview of read quality metrics, including per-base quality, GC content, and adapter contamination.",
        "Forward and reverse paired-end reads were merged using PandaSeq. The parameters used for merging were: -o', '5', '-O', '10', '-L', '192', '-l', '192', where '-o' specifies the minimum overlap length, '-O' specifies the maximum overlap length, '-l' sets the minimum read length after merging, and '-L' sets the maximum read length.",
        "The merged reads were verified for quality and integrity using FastQC. This step ensured that the merging process did not introduce significant errors or biases in the data, and that the merged reads met the quality standards required for downstream analysis.",
        "The merged reads were checked for the presence of 4 different 20-nucleotide direct repeat (DR) sequence. The presence of these DRs was verified by search their sequences in specific position in the merged reads. The total number of reads containing these DRs was reported, along with the percentage of total reads.",
        "The internal sequence with 140bp was extracted from the merged reads. This sequence is crucial for further analysis as it contains 4 guide RNAs and 3 DRs - DR4, DR3, and DR1.",
        "The internal sequence was processed to calculate the number of times the same read sequence was observed. This step is essential for identifying the most abundant sequences and understanding the diversity of the internal sequence.",
        "The internal sequence was verified against a library database to ensure that the sequences are valid and correspond to known sequences. This step is crucial for confirming the accuracy of the internal sequence and its relevance to the study.",

    ]

    for title, step in zip(pipeline_title, pipeline_steps):
        story.append(Paragraph(title, styles['h3']))
        story.append(Paragraph(step, styles['ReportBodyText']))
        story.append(Spacer(1, 0.1 * inch))

    story.append(PageBreak())  # Start detailed sections on a new page

    return story