# Bionl_Lean_Call

Germline variant calling and reporting for ACMG SF genes, including QC and functional annotation—designed for research, not diagnostic use.

---

## Inputs

### Sample Sheet (CSV)

Provide a CSV file listing your samples and their data files.

**Example:**

```csv
patient,sex,status,sample,fastq_1,fastq_2,lane
Patient001,XX,0,Sample_A,data/sample_A_R1.fastq.gz,data/sample_A_R2.fastq.gz,1
Patient002,XY,0,Sample_B,data/sample_B_R1.fastq.gz,data/sample_B_R2.fastq.gz,1
Patient003,XX,1,Sample_C,data/sample_C_R1.fastq.gz,data/sample_C_R2.fastq.gz,1
```

The CSV should include:
- **patient**: Patient identifier
- **sex**: Sex (XX for female, XY for male, or 0 if unknown)
- **status**: Affected status (0 for unaffected, 1 for affected)
- **sample**: Sample identifier (must be unique)
- **fastq_1**: Path to forward reads file
- **fastq_2**: Path to reverse reads file
- **lane**: Sequencing lane number

### Output Directory

Specify where you want results saved. The pipeline will create organized folders for each sample.

---

## Outputs

Results are organized by sample in the output directory:

```
outdir/
└── {SAMPLE}/
    ├── vcf/          # Variant files (filtered and annotated)
    ├── qc/           # Quality metrics and coverage statistics
    └── reports/      # Excel and HTML reports
```

**What you'll find:**

- **VCF files**: Identified genetic variants with annotations
- **QC metrics**: 
  - Coverage analysis (20x, 30x, 50x, 100x thresholds)
  - Coverage gaps identification
  - Alignment statistics (mapped reads, duplicates)
  - Read balance metrics (R1/R2 ratios, strand balance)
  - Sex determination from sequencing data
  - Variant statistics (transition/transversion ratios)
- **Reports**: Easy-to-read Excel summaries with variant classifications and clinical annotations

---

## How to Run

### Using the Platform UI

1. Navigate to the pipeline in your platform interface
2. Upload or select your **sample sheet** (CSV file)
3. Specify your **output directory** path
4. Click **Run**

### Command Line Alternative

```bash
nextflow run main.nf \
  --samplesheet samples.csv \
  --outdir results
```

That's it! The pipeline handles everything else automatically.

---

## Support

For questions or assistance:
- Contact your bioinformatics team
- Email: support@example.com

**Pipeline version**: 1.0.0
