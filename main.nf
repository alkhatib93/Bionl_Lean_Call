nextflow.enable.dsl=2

// ── Import Sarek pieces from the vendored repo (per your grep paths) ───────────
include { PIPELINE_INITIALISATION; PIPELINE_COMPLETION } \
  from './external/sarek/subworkflows/local/utils_nfcore_sarek_pipeline'

include { NFCORE_SAREK } \
  from './external/sarek/main.nf'

// ── (Optional) Import your post-Sarek mini-pipeline (e.g. VEP) ─────────────────
include { POST_SAREK } from './modules/vep.nf'

// ── Params (defaults). Map your samplesheet → Sarek's expected --input ─────────
//params.samplesheet   = params.samplesheet   ?: "${workflow.projectDir}/data/samplesheet.csv"
params.input = params.input ?: params.samplesheet
//params.remove('samplesheet')
params.outdir = params.outdir ?: params.output
//params.outdir  = params.outdir  ?: "${workflow.projectDir}/results/sarek"
params.bed     = params.bed     ?: "${workflow.projectDir}/data/annotated_merged_MANE_deduped.bed"
//params.intervals     = params.intervals     ?: "${workflow.projectDir}/data/chr22_targets.bed"
params.run_sarek = (params.run_sarek instanceof Boolean) ? params.run_sarek : true
// Any extra Sarek params you'll pass on CLI (e.g., --genome, --fasta, --tools, …)
// ── Fail fast if missing ──
if( params.run_sarek ) {
  if( !params.input )  error "Missing --input (samplesheet CSV) when run_sarek=true"
  if( !params.outdir ) error "Missing --outdir when run_sarek=true"
} else {
  //if( !params.sarek_outdir ) error "Missing --sarek_outdir when run_sarek=false"
  if( !params.post_samplesheet && !params.sarek_outdir )
    error "When run_sarek=false provide either --post_samplesheet or --sarek_outdir"
  if( params.post_samplesheet && params.sarek_outdir )
    error "Cannot provide both --post_samplesheet and --sarek_outdir. Choose one."
}
// ── Helper workflow to collect Sarek outputs ────────────────────────────────────
workflow COLLECT_SAREK_OUTPUTS {
  take:
    trigger    // A channel that acts as a completion signal
    outdir     // The output directory to search

  main:
    def isGCS = outdir.toString().startsWith('gs://')
    // Use trigger to wait, then collect files
    vcf_ch = trigger
      .flatMap { 
        file("${outdir}/variant_calling/*/*/*.vcf.gz", checkIfExists: isGCS)
      }
      .filter { vcf -> 
        vcf.name.endsWith('.vcf.gz') && 
        !vcf.name.contains('.g.vcf.gz') && 
        !vcf.name.endsWith('.tbi') 
      }
      .map { vcf -> tuple(vcf.parent.name, vcf) }

    bam_ch = trigger
      .flatMap { 
        file("${outdir}/preprocessing/mapped/*/*.sorted.bam", checkIfExists: isGCS)
      }
      .map { bam -> 
        def sample = bam.parent.name
        def bai = file("${bam}.bai")
        if (!bai.exists()) {
          bai = file("${bam.parent}/${bam.baseName}.bai")
        }
        if (!bai.exists()) {
          error "BAM index not found for ${bam}"
        }
        tuple(sample, bam, bai)
      }

  emit:
    vcf = vcf_ch
    bam = bam_ch
}

// ── Top-level pipeline ─────────────────────────────────────────────────────────
workflow {

  // ── BED channel (needed for all scenarios) ──
  def bedFile = params.bed ? file(params.bed) : null
  if (!bedFile?.exists()) error "BED file not found: ${params.bed}"
  bed_ch = Channel.value(bedFile)

  if (params.sarek_outdir) {
    log.info ">>> Skipping Sarek run, using existing results in ${params.sarek_outdir}"
    log.info ">>> Checking if ${params.sarek_outdir} is a GCS bucket"
    def isGCS = params.sarek_outdir.toString().startsWith('gs://')
    // ── Collect VCFs ──
    vcf_ch = Channel
      .fromPath("${params.sarek_outdir}/variant_calling/*/*/*.vcf.gz", checkIfExists: isGCS)
      .filter { vcf -> 
        vcf.name.endsWith('.vcf.gz') && 
        !vcf.name.contains('.g.vcf.gz') && 
        !vcf.name.endsWith('.tbi') 
      }
      .map { vcf -> 
        def sample = vcf.parent.name
        tuple(sample, vcf) 
      }

    // ── Collect BAMs with BAI ──
    bam_ch = Channel
      .fromPath("${params.sarek_outdir}/preprocessing/mapped/*/*.sorted.bam", checkIfExists: isGCS)
      .map { bam -> 
        def sample = bam.parent.name
        def bai = file("${bam}.bai")
        if (!bai.exists()) {
          bai = file("${bam.parent}/${bam.baseName}.bai")
        }
        if (!bai.exists()) {
          error "BAM index not found for ${bam}. Expected ${bam}.bai or ${bam.baseName}.bai"
        }
        tuple(sample, bam, bai) 
      }
    
    vcf_ch.view { s, v -> "VCF -> ${s} :: ${v}" }
    bam_ch.view { s, a, i -> "ALN -> ${s} :: ${a} | IDX ${i}" }
    // Optional safety checks
    vcf_ch.ifEmpty { error "No VCFs found in ${params.sarek_outdir}/variant_calling/*/*/*.vcf.gz" }
    bam_ch.ifEmpty  { error "No BAMs found in ${params.sarek_outdir}/preprocessing/mapped/*/*.sorted.bam" }

    // ── Run post-Sarek steps ──
    POST_SAREK(vcf_ch, bam_ch, bed_ch)

  } else if (params.post_samplesheet) {
    log.info ">>> Running post-Sarek from custom samplesheet ${params.post_samplesheet}"

    // Read samplesheet and validate files exist
    Channel
      .fromPath(params.post_samplesheet, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        def v = file(row.vcf)
        def b = file(row.bam)
        def bi = file(row.bai ?: "${b}.bai")
        if (!v.exists()) error "VCF not found: ${v}"
        if (!b.exists()) error "BAM not found: ${b}"
        if (!bi.exists()) {
          bi = file("${b.parent}/${b.baseName}.bai")
          if (!bi.exists()) error "BAI not found for ${b}"
        }
        tuple(row.sample, v, b, bi)
      }
      .multiMap { sample, vcf, bam, bai ->
        vcf: tuple(sample, vcf)
        bam: tuple(sample, bam, bai)
      }
      .set { result }

    vcf_ch = result.vcf
    bam_ch = result.bam

    // ── Run post-Sarek steps ──
    POST_SAREK(vcf_ch, bam_ch, bed_ch)

  } else {
    log.info ">>> Running Sarek as part of the pipeline"

    PIPELINE_INITIALISATION(
      params.version,
      params.validate_params,
      params.monochrome_logs,
      args,
      params.outdir,
      params.input
    )

    NFCORE_SAREK(PIPELINE_INITIALISATION.out.samplesheet)

    PIPELINE_COMPLETION(
      params.email,
      params.email_on_fail,
      params.plaintext_email,
      params.outdir,
      params.monochrome_logs,
      params.hook_url,
      NFCORE_SAREK.out.multiqc_report
    )

    // Collect outputs after Sarek completes
    COLLECT_SAREK_OUTPUTS(
      NFCORE_SAREK.out.multiqc_report,
      params.outdir
    )

    // ── Run post-Sarek steps ──
    POST_SAREK(
      COLLECT_SAREK_OUTPUTS.out.vcf, 
      COLLECT_SAREK_OUTPUTS.out.bam, 
      bed_ch
    )
  }
}