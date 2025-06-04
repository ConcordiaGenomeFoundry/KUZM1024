from reportlab.platypus import Paragraph, Spacer, PageBreak
from reportlab.lib.units import inch


def generate_report(story, styles):
    story.append(Paragraph("<b>Sequencing Data Pipeline</b>", styles['ReportSectionHeading']))
    story.append(Paragraph(
        "The following steps outline the workflow to process raw illumina sequencing reads through to alignment, filtering and data analysis used to generate this report:",
        styles['ReportBodyText']))
    story.append(Spacer(1, 0.3 * inch))

    pipeline_title = [
        "1. Read Quality Assessment",
        "2. Read Merging",
        "3. FASTQ to FASTA Conversion",
        "4. Database Preparation",
        "5. Sequence Alignment",
        "6. Alignment Filtering"
    ]

    pipeline_steps = [
        "Initial quality control of raw sequencing reads was performed using FastQC to generate individual read quality reports. These reports were then aggregated and summarized with MultiQC for a comprehensive overview of read quality metrics, including per-base quality, GC content, and adapter contamination.",
        "Forward and reverse paired-end reads were merged using PEAR (Paired-End reAd merger). The merging process was executed with a minimum quality score threshold of 25 (-q 25), a minimum overlap of 10 bases and minimun length od reads of 20 (-t 20).",
        "Merged reads, originally in FASTQ format, were converted to FASTA format. This step prepares the sequences for database creation and subsequent alignment.",
        "A custom sequence database was created in FASTA format. This database serves as the reference for the subsequent alignment step.",
        "Reads were aligned against the prepared database using MMseqs2. This tool was chosen for its efficiency in performing sensitive sequence similarity searches.",
        "A filtering step was applied to the alignment results. Only alignments meeting specific criteria were retained: "
            "<li><strong>Identity</strong>: Alignments with less than 97% (3 mismatches) sequence identity were removed.</li>"
            "<li><strong>Length</strong>: Alignments not precisely 140 base pairs (bp) in length were excluded to ensure consistent sequence representation with the gRNAs length.</li>"
    ]

    for title, step in zip(pipeline_title, pipeline_steps):
        story.append(Paragraph(title, styles['h3']))
        story.append(Paragraph(step, styles['ReportBodyText']))
        story.append(Spacer(1, 0.1 * inch))

    story.append(PageBreak())  # Start detailed sections on a new page

    return story