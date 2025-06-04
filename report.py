from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.units import inch
from reportlab.lib import colors
import numpy as np
import pandas as pd
import glob, os, io
import matplotlib
matplotlib.use('Agg') # This line must be before importing pyplot
import matplotlib.pyplot as plt
import main  # Assuming 'main' module is in the same directory and accessible
import pipeline_info


def calculate_gini(grna_counts):
    """
    :param grna_counts: List or array of gRNA counts.
    :return: Gini index as a float.
    Calculates the Gini index for a given set of gRNA counts.
    """

    grna_array = np.array(grna_counts)

    """Compute Gini coefficient of array of values"""
    diffsum = 0
    for i, count in enumerate(grna_array[:-1], 1):
        # For each count, calculate the absolute differences with all subsequent counts
        diffsum += np.sum(np.abs(count - grna_array[i:]))
    # Normalize the Gini index
    gini = diffsum / (len(grna_array) ** 2 * np.mean(grna_array))
    return gini


def create_gini_chart(gini_index):
    """Generates a bar chart visualizing the Gini index."""
    fig, ax = plt.subplots(figsize=(6, 3))  # Adjusted for narrower space

    # Create a bar for the Gini index value
    ax.bar(['Gini Index'], [gini_index], color='#336699')  # Corrected: Pass hex string directly

    # Add a horizontal line at 0.5 for reference
    ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8, label='Midpoint (0.5)')

    # Add text label for the Gini index value on the bar
    ax.text(0, gini_index + 0.02, f'{gini_index:.4f}', ha='center', va='bottom', fontsize=10, color='black')

    ax.set_ylim(0, 1)  # Gini index is between 0 and 1
    ax.set_ylabel('Value')
    ax.set_title('Gini Index (Uniformity)', fontsize=12)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(loc='upper right', frameon=False, fontsize=8)  # Add legend for the midpoint line

    # Save to a BytesIO object (in memory)
    buf = io.BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    plt.close(fig)  # Close the plot to free memory
    return buf


def create_coverage_chart(coverage):
    """Generates a simple bar chart for coverage."""
    fig, ax = plt.subplots(figsize=(6, 3))

    ax.bar(['Coverage'], [coverage], color='#28a745')  # Corrected: Pass hex string directly
    ax.set_ylim(0, 100)
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Library Coverage', fontsize=12)
    ax.text(0, coverage + 1, f'{coverage:.2f}%', ha='center', va='bottom', fontsize=10, color='black')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    buf = io.BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    plt.close(fig)
    return buf


def create_percentage_reads_chart(percentage_at_least_5_reads):
    """Generates a pie chart for % gRNAs with >=5 reads."""
    fig, ax = plt.subplots(figsize=(5, 3))  # Slightly taller for pie chart

    labels = ['≥5 Reads', '<5 Reads']
    sizes = [percentage_at_least_5_reads, 100 - percentage_at_least_5_reads]
    colors_pie = ['#007bff', '#6c757d']  # Corrected: Pass hex strings directly

    # Only plot if sizes are valid (e.g., if total counts are zero, sizes might be [0, 0])
    if sum(sizes) > 0:
        ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors_pie,
               wedgeprops={'edgecolor': 'white'})
    else:
        # Handle case where no reads are detected, or counts are 0
        ax.text(0, 0, 'No gRNAs detected', ha='center', va='center', fontsize=12, color='gray')
        ax.set_xticks([])
        ax.set_yticks([])

    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax.set_title('% gRNAs with ≥5 Reads', fontsize=12)

    buf = io.BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    plt.close(fig)
    return buf


def create_single_boxplot_with_mean_std(data, dataset_label="gRNA Read Counts"):
    """
    Creates a single boxplot with the mean and standard deviation indicated,
    returning it as a BytesIO buffer and a list of outliers.

    Args:
        data (array-like): The dataset for which to create the boxplot.
        dataset_label (str): The label for the dataset on the x-axis.
    """
    if not isinstance(data, (list, np.ndarray)):
        data = np.array(data)

    if data.size == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, 'No data available for Boxplot', ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f'Boxplot: {dataset_label}')
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=300)
        buf.seek(0)
        plt.close(fig)
        return buf, [] # Return empty list for outliers if no data

    data_mean = np.mean(data)
    data_std = np.std(data)

    fig, ax = plt.subplots(figsize=(6, 4))

    boxplot_elements = ax.boxplot(data, vert=True, patch_artist=True, widths=0.3, showfliers=True)

    for patch in boxplot_elements['boxes']:
        patch.set_facecolor('#ADD8E6')
        patch.set_edgecolor('black')

    for median in boxplot_elements['medians']:
        median.set(color='red', linewidth=2)

    for whisker in boxplot_elements['whiskers']:
        whisker.set(color='gray', linestyle='-')

    for cap in boxplot_elements['caps']:
        cap.set(color='gray', linewidth=2)

    for flier in boxplot_elements['fliers']:
        flier.set(marker='o', color='gray', alpha=0.5)

    ax.set_yscale('log')

    ax.plot(1, data_mean, marker='D', markersize=5, color='#DC143C', label='Mean')
    ax.vlines(1, data_mean - data_std, data_mean + data_std, color='#DC143C', linestyle='-', linewidth=2, label='Std Dev')

    ax.set_xticks([1])
    ax.set_xticklabels([dataset_label])
    ax.set_ylabel('Count')
    ax.set_title('gRNA Read Count Distribution', fontsize=12)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    plt.close(fig)

    # Extract outliers
    outliers = []
    # The 'fliers' element in boxplot_elements is a list of Line2D artists.
    # Each artist's ydata contains the outlier values.
    for flier in boxplot_elements['fliers']:
        outliers.extend(flier.get_ydata())

    return buf, outliers


class ReportlabHelper:
    """Helper class for ReportLab elements."""

    def __init__(self):
        self.styles = getSampleStyleSheet()
        self.styles.add(ParagraphStyle(name='ReportTitle',
                                       parent=self.styles['Title'],
                                       fontSize=24,
                                       spaceAfter=14,
                                       alignment=TA_CENTER,
                                       textColor=colors.HexColor('#003366')))  # Dark blue for title
        self.styles.add(ParagraphStyle(name='ReportSubTitle',
                                       parent=self.styles['Title'],
                                       fontSize=18,
                                       spaceAfter=10,
                                       alignment=TA_CENTER,
                                       textColor=colors.HexColor('#003366')))  # Dark blue for title
        self.styles.add(ParagraphStyle(name='ReportSectionHeading',
                                       parent=self.styles['h1'],
                                       fontSize=16,
                                       spaceBefore=12,
                                       spaceAfter=6,
                                       textColor=colors.HexColor('#336699')))  # Medium blue for headings
        self.styles.add(ParagraphStyle(name='ReportBodyText',
                                       parent=self.styles['Normal'],
                                       fontSize=10,
                                       leading=14,
                                       spaceAfter=6))
        self.styles.add(ParagraphStyle(name='FooterStyle',
                                       parent=self.styles['Normal'],
                                       fontSize=8,
                                       alignment=TA_CENTER,
                                       textColor=colors.HexColor('#666666')))

    def get_styles(self):
        return self.styles

    def header_footer(self, canvas, doc):
        """Adds a header and footer to each page."""
        canvas.saveState()
        # Header
        logo_path = "logo.svg"  # Make sure you have a logo.png in your directory
        if os.path.exists(logo_path):
            logo = Image(logo_path, width=0.75 * inch, height=0.75 * inch)
            logo.wrapOn(canvas, doc.width, doc.topMargin)
            logo.drawOn(canvas, doc.leftMargin, doc.height + doc.topMargin - 0.8 * inch)

        # Footer
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(colors.HexColor('#666666'))
        canvas.drawCentredString(doc.pagesize[0] / 2.0, 0.75 * inch, f"Page {doc.page}")
        canvas.drawString(doc.leftMargin, 0.75 * inch, "Generated: " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"))
        canvas.restoreState()


def generate_report(path_sequencing_data, merged_reads_folder, counts_txt_file):
    """
    Generates a PDF report summarizing sequencing analysis metrics.

    :param path_sequencing_data: Path to the sequencing data directory.
    :param merged_reads_folder: Name of the folder containing merged reads.
    :return: A report in pdf summarizing the coverage, uniformity, fold-80, mean reads, and % gRNAs.
    """
    output_filename = "sequencing_analysis_report.pdf"

    try:
        df_grna_count = pd.read_csv(counts_txt_file, sep="\t", header=None, names=["qseqid", "count"])
        grna_counts = df_grna_count["count"].values
    except Exception as e:
        print(f"Error reading {counts_txt_file}: {e}")
        return

    fasta_dir = os.path.join(path_sequencing_data, merged_reads_folder, 'fasta')
    fasta_files = glob.glob(os.path.join(fasta_dir, '*.fasta'))
    fasta_file_base = "N/A"
    if fasta_files:
        fasta_file_base = os.path.basename(fasta_files[0]).split('.')[0]

    # Calculate metrics for the summary table
    total_grnas, average_length = 0, 0
    try:
        total_grnas, average_length = main.count_reads(path_sequencing_data)
    except Exception as e:
        print(f"Warning: Could not get total gRNAs from main.count_reads: {e}")

    detected_grnas = np.sum(grna_counts > 0)
    coverage = (detected_grnas / total_grnas) * 100 if total_grnas > 0 else 0
    gini_index = calculate_gini(grna_counts)
    mean_reads = np.mean(grna_counts)
    grnas_at_least_5_reads = np.sum(grna_counts >= 5)
    percentage_at_least_5_reads = (grnas_at_least_5_reads / len(grna_counts)) * 100 if len(grna_counts) > 0 else 0

    # ---
    # Charts
    boxplot_buffer, outliers = create_single_boxplot_with_mean_std(grna_counts, dataset_label="gRNA Read Counts")
    completeness_chart_buffer = create_coverage_chart(coverage)
    gini_chart_buffer = create_gini_chart(gini_index)
    percentage_reads_chart_buffer = create_percentage_reads_chart(percentage_at_least_5_reads)

    # Initializing Reportlab document
    doc = SimpleDocTemplate(output_filename, pagesize=letter,
                            leftMargin=0.8 * inch,
                            rightMargin=0.8 * inch,
                            topMargin=1.0 * inch,
                            bottomMargin=0.75 * inch)

    helper = ReportlabHelper()
    styles = helper.get_styles()
    story = []

    # Logo
    logo_path = "logo.png"  # Ensure you have a logo.svg in your directory
    if os.path.exists(logo_path):
        # reduce the image but keep aspect ratio
        logo = Image(logo_path, width=4.2 * inch, height=1.5 * inch)
        logo.hAlign = 'CENTER'
        story.append(logo)

    story.append(Spacer(1, 2 * inch))  # Vertical space
    story.append(Paragraph(f"Sequencing Analysis Report", styles['ReportTitle']))
    story.append(Paragraph(f"for {fasta_file_base}", styles['ReportSubTitle']))  # Subtitle for specific file
    story.append(Spacer(1, 0.5 * inch))
    story.append(Paragraph("Concordia Genome Foundry", styles['h3']))
    story.append(Paragraph("Date: " + pd.Timestamp.now().strftime("%B %d, %Y"), styles['h3']))
    story.append(PageBreak())

    # Pipeline information
    story = pipeline_info.generate_report(story, styles)

    # Main Content Pages
    # The header and footer will be applied to these pages automatically

    # ---
    # 1. Summary Table
    story.append(Paragraph("<b>Sequencing Data Analysis</b>", styles['ReportSectionHeading']))
    story.append(Paragraph("<b>Summary</b>", styles['Heading2']))
    story.append(Spacer(1, 0.1 * inch))
    summary_data = [
        ["Metric", "Value"],
        ["Total gRNA constructs in library (from template)", f"{total_grnas:,}"],
        ["Detected gRNA constructs (>0 reads)", f"{detected_grnas:,}"],
        ["Total of outliers detected", f"{len(outliers):,}"],
        ["Missing gRNA constructs", f"{(total_grnas - detected_grnas):,}"],
        ["Library coverage", f"{coverage:.2f}%"],
        ["Gini Index (Uniformity)", f"{gini_index:.4f}"],
        ["Mean reads per gRNA", f"{mean_reads:,.2f}"],
        ["Total of gRNAs with ≥5 Reads", f"{grnas_at_least_5_reads:,}"],
        ["% gRNAs with ≥5 Reads", f"{percentage_at_least_5_reads:.2f}%"],
    ]

    summary_table = Table(summary_data, colWidths=[3.5 * inch, 2.5 * inch])
    summary_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E6EEF6')),  # Header background
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.HexColor('#003366')),  # Header text color
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#DDDDDD')),
        ('BOX', (0, 0), (-1, -1), 1, colors.HexColor('#AAAAAA')),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('TOPPADDING', (0, 1), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
    ]))
    story.append(summary_table)
    story.append(Spacer(1, 0.4 * inch))

    story.append(Paragraph(
        "This report provides a comprehensive analysis of the sequencing data, focusing on key metrics to assess "
        "library quality and coverage uniformity. File used to calculate these metrics was <strong>"
        + str(counts_txt_file) + ".</strong> The following sections detail each metric.",
        styles['BodyText']))
    story.append(PageBreak())  # Start detailed sections on a new page

    # ---
    # 2. Detailed Analysis Sections

    # gRNA Read Count Distribution Boxplot
    story.append(Paragraph("<b>1. gRNA Read Count Distribution</b>", styles['ReportSectionHeading']))
    story.append(Paragraph(
        "This boxplot visualizes the distribution of read counts per gRNA, showing the median, interquartile range "
        "(IQR), and potential outliers. The red dot indicates the mean, and the red vertical line represents one "
        "standard deviation above and below the mean. List of outliers: <strong> outliers.csv. </strong>",
        styles['ReportBodyText']))

    story.append(Image(boxplot_buffer, width=5.0 * inch, height=4.0 * inch))  # Adjust size as needed
    story.append(Spacer(1, 0.2 * inch))

    # If outliers were detected, save them to a CSV file
    if outliers:
        # Create a csv file with only the outliers using the df_grna_count file
        outliers_df = df_grna_count[df_grna_count['count'].isin(outliers)]
        outliers_csv_path = "outliers.csv"
        outliers_df.to_csv(outliers_csv_path, index=False)

    # Library coverage
    story.append(Paragraph("<b>2. Library coverage</b>", styles['ReportSectionHeading']))
    story.append(Paragraph(
        f"Library coverage measures the proportion of expected gRNAs that were detected in the sequencing run. "
        f"A higher Library coverage percentage indicates a more comprehensive representation of the gRNA library.",
        styles['BodyText']))
    story.append(Paragraph(f"Total gRNAs (from template): <b>{total_grnas:,}</b>", styles['BodyText']))
    story.append(Paragraph(f"Detected gRNAs (with >0 reads): <b>{detected_grnas:,}</b>", styles['BodyText']))
    story.append(Paragraph(f"<b>Library coverage: {coverage:.2f}%</b>", styles['BodyText']))

    story.append(Image(completeness_chart_buffer, width=4.5 * inch, height=2.5 * inch))  # Adjust size as needed
    story.append(Spacer(1, 0.2 * inch))

    # Uniformity (Gini Index)
    story.append(Paragraph("<b>3. Uniformity (Gini Index)</b>", styles['ReportSectionHeading']))
    story.append(Paragraph(
        f"The Gini index quantifies the uniformity of gRNA distribution. A value closer to 0 indicates perfect "
        f"uniformity (all gRNAs have equal reads), while a value closer to 1 indicates high inequality "
        f"(a few gRNAs dominate the reads). For library screening, a lower Gini index is generally desired.",
        styles['BodyText']))
    story.append(Paragraph(f"<b>Gini Index: {gini_index:.4f}</b>", styles['BodyText']))

    # Add Gini Index chart
    story.append(Image(gini_chart_buffer, width=4.5 * inch, height=2.5 * inch))  # Adjust size as needed
    story.append(Spacer(1, 0.2 * inch))

    # % gRNAs with ≥5 Reads
    story.append(Paragraph("<b>4. % gRNAs with ≥5 Reads</b>", styles['ReportSectionHeading']))
    story.append(Paragraph(
        f"This metric indicates the percentage of unique gRNAs that have received at least 5 sequencing reads. "
        f"This threshold is often used to assess the effective coverage of the library, as very low read counts "
        f"might not be sufficient for downstream analysis.",
        styles['BodyText']))
    story.append(Paragraph(f"gRNAs with ≥5 reads: <b>{grnas_at_least_5_reads:,}</b>", styles['BodyText']))
    story.append(
        Paragraph(f"<b>Percentage of gRNAs with ≥5 reads: {percentage_at_least_5_reads:.2f}%</b>", styles['BodyText']))

    story.append(Image(percentage_reads_chart_buffer, width=4 * inch, height=2.5 * inch))
    story.append(Spacer(1, 0.2 * inch))

    # Top 20 gRNAs Table (NEW SECTION)
    story.append(Paragraph("<b>5. Top 20 gRNAs by Read Count</b>", styles['ReportSectionHeading']))
    story.append(Paragraph(
        "This table lists the top 20 gRNAs with the highest number of sequencing reads. This highlights potential "
        "biases or highly represented gRNAs within the library.",
        styles['BodyText']))

    # Prepare data for the top 20 gRNAs table
    if not df_grna_count.empty:
        top_20_grnas = df_grna_count.sort_values(by="count", ascending=False).head(20)
        top_20_data = [["Rank", "gRNA ID", "Read Count"]]
        for i, (index, row) in enumerate(top_20_grnas.iterrows()):
            top_20_data.append([i + 1, row["qseqid"], f"{row['count']:,}"])

        top_20_table = Table(top_20_data, colWidths=[0.8 * inch, 3.5 * inch, 1.7 * inch])
        top_20_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E6EEF6')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.HexColor('#003366')),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#DDDDDD')),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor('#AAAAAA')),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('TOPPADDING', (0, 1), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
        ]))
        story.append(top_20_table)
    else:
        story.append(Paragraph("No gRNA read count data available to generate top 20 list.", styles['BodyText']))
    story.append(Spacer(1, 0.2 * inch))

    # ---
    # Footer
    story.append(Spacer(1, 0.4 * inch))
    story.append(Paragraph(
        "This report was created for C_KUZM1024 using a Illumina sequencing analysis pipeline. For any questions or "
        "further analysis, please contact Flavia Araujo or Jing Cheng. Source code available "
        "at: https://github.com/ConcordiaGenomeFoundry/KUZM1024",
        styles['FooterStyle']))

    # ---
    # Build the PDF
    try:
        doc.build(story)
        print(f"Report '{output_filename}' created successfully!")
    except Exception as e:
        print(f"Error creating PDF: {e}")