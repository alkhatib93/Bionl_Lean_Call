# Bionl_Lean_Call

Germline variant calling and reporting for ACMG SF genes, including QC and functional annotation—designed for research, not diagnostic use.

---

## Inputs

The pipeline supports two input modes:

### Option 1: Run from FASTQ files (Full Pipeline)

Provide a CSV file listing your samples and their FASTQ files.

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

**OR**

### Option 2: Run from Sarek Output (Post-Processing Only)

If you already have Sarek results (BAM and VCF files), use a post-processing sample sheet:

**Example (`post_samplesheet.csv`):**

```csv
sample,vcf,bam,bai
HG003,/path/to/HG003.deepvariant.vcf.gz,/path/to/HG003.sorted.bam,/path/to/HG003.sorted.bam.bai
```

The CSV should include:
- **sample**: Sample identifier
- **vcf**: Path to the VCF file (can be from DeepVariant, FreeBayes, etc.)
- **bam**: Path to the aligned BAM file
- **bai**: Path to the BAM index file

You can provide either:
- `--sarek_outdir`: Path to a Sarek output directory (automatically finds VCF/BAM files)
- `--post_samplesheet`: Path to a CSV with explicit file paths (use this for more control)

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
2. Choose your input mode:
   - **Full pipeline**: Upload or select your **FASTQ sample sheet** (CSV file)
   - **Post-processing**: Upload **post_samplesheet.csv** OR specify **Sarek output directory**
3. Specify your **output directory** path
4. Click **Run**

### Command Line Examples

**Full pipeline (from FASTQ files):**
```bash
nextflow run main.nf \
  --samplesheet samples.csv \
  --outdir results
```

**Post-processing from Sarek output directory:**
```bash
nextflow run main.nf \
  --sarek_outdir /path/to/sarek/results \
  --outdir results
```

**Post-processing with explicit file paths:**
```bash
nextflow run main.nf \
  --post_samplesheet data/post_samplesheet.csv \
  --outdir results
```

That's it! The pipeline handles everything else automatically.

---

## Support

For questions or assistance:
- Contact your bioinformatics team
- Email: support@example.com

**Pipeline version**: 1.0.0
