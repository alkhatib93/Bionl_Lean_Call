# Usage Guide

## Overview

The Bionl_Lean_Call analyzes genomic sequencing data to identify variants and assess data quality. Simply provide a sample sheet listing your samples and an output directory location. The pipeline processes each sample independently and generates comprehensive reports with variant annotations and quality metrics.

Results include filtered variant calls, detailed quality control statistics, and user-friendly Excel reports suitable for clinical review.

---

## Input Samplesheet

Create a CSV file with your sample information:

```csv
patient,sex,status,sample,fastq_1,fastq_2,lane
PT001,XX,0,SampleA,/data/sampleA_R1.fastq.gz,/data/sampleA_R2.fastq.gz,1
PT002,XY,0,SampleB,/data/sampleB_R1.fastq.gz,/data/sampleB_R2.fastq.gz,1
PT003,XX,1,SampleC,/data/sampleC_R1.fastq.gz,/data/sampleC_R2.fastq.gz,1
```

**Columns:**
- `patient` - Patient identifier
- `sex` - Sex (XX for female, XY for male, or 0 if unknown)
- `status` - Affected status (0 for unaffected, 1 for affected)
- `sample` - Sample identifier (must be unique)
- `fastq_1` - Full path to forward reads
- `fastq_2` - Full path to reverse reads
- `lane` - Sequencing lane number (typically 1)

---

## Outputs at a Glance

```
results/
├── SampleA/
│   ├── vcf/
│   │   └── SampleA.vep_annotated.vcf.gz
│   ├── qc/
│   │   ├── SampleA_coverage_summary.txt
│   │   ├── SampleA_flagstat.txt
│   │   └── SampleA.mosdepth.summary.txt
│   └── reports/
│       ├── SampleA_variants_lean_v1.xlsx
│       └── SampleA_clinical_report.html
├── SampleB/
│   └── ...
└── SampleC/
    └── ...
```

Each sample gets three main folders:
- **vcf/** - Variant call files with annotations
- **qc/** - Quality control metrics including:
  - Coverage statistics at multiple thresholds (10x, 20x, 30x, 50x, 100x)
  - Coverage gaps below critical thresholds
  - Alignment quality (total reads, mapped reads, duplicate rates)
  - Read balance analysis (R1/R2 ratios, forward/reverse strand balance)
  - Sex determination from X/Y chromosome coverage
  - Variant quality statistics
- **reports/** - Excel and HTML summary reports

---

## Run

### Via Platform UI

1. Open the pipeline in your platform
2. **Upload samplesheet**: Select your CSV file
3. **Set output directory**: Choose where to save results (e.g., `/projects/analysis_2024`)
4. Click **Start Pipeline**

The platform will show progress as samples are processed.

### Via Command Line

If running outside the platform, use this simple command:

```bash
nextflow run main.nf \
  --samplesheet my_samples.csv \
  --outdir my_results
```

---

## Notes

**Understanding Report Messages:**

When reviewing your Excel reports, you may see statements like:

- **"No variants detected in this category"** - This means no variants matched the specific filters for that section (e.g., no pathogenic variants found). This is often a normal result.

- **"No exons below 20x coverage"** - This is good news! It means all target regions have sufficient sequencing depth.

- **"No coverage gaps identified"** - All areas of interest were adequately sequenced.

These "no data" statements typically indicate quality results rather than errors.

**Quality Control Files:**

The QC folder contains detailed metrics files:
- **Coverage files** - Per-region coverage at different depth thresholds
- **Gap files** - Regions below 20x or 30x coverage
- **Flagstat** - Read alignment summary (total, mapped, duplicates)
- **Stats** - Detailed BAM statistics
- **R1/R2 ratios** - Read pair balance per region
- **Strand balance** - Forward/reverse strand balance per region
- **Sex check** - Predicted sex from X/Y chromosome coverage
- **Variant stats** - Transition/transversion ratios and quality distributions

**Report Sections:**

Your Excel report contains multiple tabs:
- **Sample Summary** - Overall quality metrics (coverage, alignment, ratios)
- **ACMG SF (P-LP)** - Clinically significant variants
- **Coverage gaps** - Any regions with low coverage
- **ACMG SF Genes Coverage** - Coverage statistics for key genes
- **PASS variants** - All high-quality variants

**Typical Processing Time:**

- Small cohorts (1-3 samples): 2-4 hours
- Medium cohorts (4-10 samples): 4-8 hours  
- Large cohorts (>10 samples): Plan accordingly

---

## Questions?

Contact your bioinformatics support team for assistance with results interpretation or technical issues.
