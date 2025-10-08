nextflow.enable.dsl=2
log.info ">>> Bionl_Lean_Call v1.0.3 – updated on 2025-10-07 <<<"
// Input params
params.samplesheet      = params.samplesheet      ?: "${workflow.projectDir}/data/samplesheet.csv"
params.post_samplesheet = params.post_samplesheet ?: "${workflow.projectDir}/data/post_samplesheet.csv"
params.intervals        = params.intervals        ?: "${workflow.projectDir}/data/chr22_targets.bed"
params.bed              = params.bed              ?: "data/annotated_merged_MANE_deduped.bed"
params.outdir           = params.outdir           ?: "${workflow.projectDir}/results"
params.scriptdir        = params.scriptdir        ?: "${workflow.projectDir}/scripts"
params.template_dir     = params.template_dir     ?: "${workflow.projectDir}/scripts/template-files"

// Sarek params
params.run_sarek        = params.run_sarek        ?: true
params.sarek_run_outdir  = params.sarek_run_outdir ?: "${params.outdir}/sarek"
params.sarek_outdir     = params.sarek_outdir     ?: null
params.sarek_rev        = params.sarek_rev        ?: "3.5.1"
params.sarek_profile    = params.sarek_profile    ?: "singularity"
params.sarek_config     = params.sarek_config     ?: "${workflow.projectDir}/conf/sarek_override.config"

// The extra args will usually come from nextflow.config, default to genome only
params.sarek_extra_args = params.sarek_extra_args ?: ""

// Annotation params
params.run_vep             = params.run_vep             ?: true
params.vep_fasta            = params.vep_fasta            ?: ""
params.revel_vcf            = params.revel_vcf            ?: ""
params.alpha_missense_vcf   = params.alpha_missense_vcf   ?: ""
params.clinvar_vcf          = params.clinvar_vcf          ?: ""


// --------- Helper to ensure required params ---------
def require(p, msg){ if(!p) error msg }
def SAREK_OUT_ABS = new File(params.sarek_run_outdir).getAbsolutePath()

// --------- Run Sarek as a process ---------
process RunSarek {
  container 'ghcr.io/nextflow-io/nextflow:24.10.0'
  tag "sarek"
  publishDir "${params.outdir}/sarek_logs", mode: 'copy'

  input:
  val samplesheet

  // emit a small marker, not the result dir
  output:
  path "sarek.done"

  shell:
  """
  set -euo pipefail
  mkdir -p "${SAREK_OUT_ABS}"

  extra_intervals=""
  if [[ -n "${params.intervals}" ]]; then
    extra_intervals="--intervals ${params.intervals}"
  fi

  nextflow run nf-core/sarek -r ${params.sarek_rev} \\
    -profile singularity\\
    --input ${samplesheet} \\
    --outdir ${SAREK_OUT_ABS} \\
    --genome GATK.GRCh38 \\
    --intervals ${params.intervals} \\
    --igenomes_ignore true \\
    --fasta ${params.vep_fasta} \\
    --skip_tools baserecalibrator \\
    --tools deepvariant \\
    --save_mapped true \\
    --save_output_as_bam true \\
    -c ${params.sarek_config} \\
    -resume

  touch sarek.done
  """   
}


// --------- Include your custom post-Sarek workflow ---------
include { POST_SAREK } from './modules/vep.nf'

// --------- Top-level pipeline logic ---------
workflow {
  def bed_ch = Channel.value(file(params.bed))
  
  def vcf_ch
  def bam_ch
  println 'Hello, World!'

  if (params.run_sarek) {
    // Run Sarek first
    def ss = Channel.value(file(params.samplesheet))
    def sarek_done = RunSarek(ss)
    
    // KEY FIX: Wait for Sarek completion, then collect outputs
    // Using .map to trigger collection after sarek.done exists
    vcf_ch = sarek_done
      .map { done_file ->
        // Once done, glob for VCF files
        file("${SAREK_OUT_ABS}/variant_calling/deepvariant/*/*.vcf.gz")
      }
      .flatten()
      .filter { vcf -> 
        vcf.name.endsWith('.vcf.gz') && 
        !vcf.name.contains('.g.vcf.gz') && 
        !vcf.name.endsWith('.tbi')
      }
      .map { vcf -> 
        def sample = vcf.parent.name
        tuple(sample, vcf)
      }

    bam_ch = sarek_done
      .map { done_file ->
        // Collect BAM files with their indices
        def bam_files = file("${SAREK_OUT_ABS}/preprocessing/mapped/*/*.sorted.bam")
        bam_files.collect { bam ->
          def sample = bam.parent.name
          def bai = file("${bam}.bai")
          // Also check for .bai without .bam prefix
          if (!bai.exists()) {
            bai = file(bam.toString().replaceAll(/\.bam$/, '.bai'))
          }
          if (!bai.exists()) {
            error "BAM index not found for ${bam}"
          }
          tuple(sample, bam, bai)
        }
      }
      .flatten()
      .collate(3) // Group back into tuples of 3
      .map { it -> tuple(it[0], it[1], it[2]) }

  //} else {
  //  // Use pre-existing Sarek output
  //  if (!params.sarek_outdir) {
  //    error "Must provide --sarek_outdir when --run_sarek is false"
  //  }
//
  //  vcf_ch = Channel
  //    .fromPath("${params.sarek_outdir}/variant_calling/deepvariant/*/*.vcf.gz", checkIfExists: true)
  //    .filter { vcf -> 
  //      vcf.name.endsWith('.vcf.gz') && 
  //      !vcf.name.contains('.g.vcf.gz') && 
  //      !vcf.name.endsWith('.tbi')
  //    }
  //    .map { vcf -> tuple(vcf.parent.name, vcf) }
//
  //  bam_ch = Channel
  //    .fromPath("${params.sarek_outdir}/preprocessing/mapped/*/*.sorted.bam", checkIfExists: true)
  //    .map { bam ->
  //      def sample = bam.parent.name
  //      def bai = file("${bam}.bai")
  //      if (!bai.exists()) {
  //        bai = file(bam.toString().replaceAll(/\.bam$/, '.bai'))
  //      }
  //      if (!bai.exists()) {
  //        error "BAM index not found for ${bam}"
  //      }
  //      tuple(sample, bam, bai)
  //    }
  //}
    } else {
    /*
     * POST-ONLY MODE (run_sarek = false)
     * Option A) --sarek_outdir points to a finished Sarek run
     * Option B) a CSV samplesheet with columns: sample, vcf, bam[, bai]
     */

    if ( params.sarek_outdir ) {
      // --- A) Reuse an existing Sarek results directory
      vcf_ch = Channel
        .fromPath("${params.sarek_outdir}/variant_calling/*/*/*.vcf.gz", checkIfExists: true)
        .filter { v -> v.name.endsWith('.vcf.gz') && !v.name.contains('.g.vcf.gz') && !v.name.endsWith('.tbi') }
        .map { vcf -> tuple(vcf.parent.name, vcf) }

      bam_ch = Channel
        .fromFilePairs("${params.sarek_outdir}/preprocessing/mapped/*/*.sorted.{bam,bai}", size: 2, flat: true, checkIfExists: true)
        .map { sample, files ->
          def bam = files.find { it.name.endsWith('.bam') }
          def bai = files.find { it.name.endsWith('.bai') }
          tuple(sample, bam, bai)
        }

    } else if ( params.post_samplesheet ) {
      // --- B) Use a post-only post_samplesheet with VCF/BAM paths
      // CSV columns required: sample, vcf, bam  (optional: bai)
      def rows = Channel
        .fromPath(params.post_samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
          // basic validation
          if( !row.sample || !row.vcf || !row.bam )
            error "Post-only post_samplesheet must have columns: sample,vcf,bam (optional bai)"
          // coerce to files
          def v = file(row.vcf)
          def b = file(row.bam)
          def i = row.bai ? file(row.bai) : file("${row.bam}.bai")
          if( !v.exists() ) error "VCF not found for sample ${row.sample}: ${v}"
          if( !b.exists() ) error "BAM not found for sample ${row.sample}: ${b}"
          if( !i.exists() ) error "BAI not found for sample ${row.sample}: ${i}"
          // emit a composite map; we’ll split into two channels below
          [ sample: row.sample as String, vcf: v, bam: b, bai: i ]
        }

      vcf_ch = rows.map { r -> tuple(r.sample, r.vcf) }
      bam_ch = rows.map { r -> tuple(r.sample, r.bam, r.bai) }

    } else {
      error "When --run_sarek false, provide either --sarek_outdir OR a post-only --post_samplesheet with columns sample,vcf,bam[,bai]."
    }
  }

  // Optional: Debug channels
  // vcf_ch.view { sample, vcf -> "VCF: ${sample} -> ${vcf.name}" }
  // bam_ch.view { sample, bam, bai -> "BAM: ${sample} -> ${bam.name}" }

  // Run post-Sarek workflow
  POST_SAREK(vcf_ch, bam_ch, bed_ch)
}
