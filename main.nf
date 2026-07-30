#!/usr/bin/env nextflow
/*
 * ==========================================================================
 *  5-Stage WGS Variant & Copy-Number Nextflow Pipeline
 * ==========================================================================
 *  Stage 1: Read Cleaning & Fragment Merging          (fastp)
 *  Stage 2: Dual-Mode Genome Mapping                  (bwa-meme)
 *  Stage 3: PCR Deduplication & Multi-Dimensional QC  (Picard / Qualimap)
 *  Stage 4: Joint Cohort Variant Genotyping & Filter   (GATK4)
 *  Stage 5: Control-Free Copy Number Profiling        (CNVkit)
 * ==========================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameters (override on the command line with --paramName value)
// ---------------------------------------------------------------------------
params.reads       = "raw_fastq"                              // dir containing *_R1_001.fastq / *_R2_001.fastq
params.ref         = "refs/ctropicalis_ref.fna"                // reference FASTA (bwa-meme indexed)
params.gtf         = "refs_gtf/GCF_000006335.3_ASM633v3_genomic.gff"
params.intervals   = "refs/genomic_intervals.list"
params.exclude     = "mk22"                                    // failed control batch, excluded from Stage 3 onward
params.outdir      = "."
params.cnvkit_exec = "cnvkit.py"

// ---------------------------------------------------------------------------
// Reference / support files
// ---------------------------------------------------------------------------
ref_fasta     = file(params.ref)
ref_index_ch  = Channel.fromPath("${params.ref}.*").collect()  // bwa-meme index files
gtf_file      = file(params.gtf)
intervals_file = file(params.intervals)

log.info """
 5-STAGE WGS PIPELINE
 ---------------------------
 reads       : ${params.reads}
 reference   : ${params.ref}
 gtf         : ${params.gtf}
 intervals   : ${params.intervals}
 exclude     : ${params.exclude}
 outdir      : ${params.outdir}
 """

// ---------------------------------------------------------------------------
// STAGE 1: Read Cleaning & Fragment Merging (fastp)
// ---------------------------------------------------------------------------
process FASTP {
    tag "${sample}"
    publishDir "${params.outdir}/01_fastp_clean", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_merged.fastq"), path("${sample}_unmerged_R1.fastq"), path("${sample}_unmerged_R2.fastq"), emit: all_reads
    path "${sample}_singleton_R1.fastq"
    path "${sample}_singleton_R2.fastq"
    path "${sample}_fastp.html"
    path "${sample}_fastp.json"

    script:
    """
    export LC_ALL=C
    fastp \
      -i "${r1}" \
      -I "${r2}" \
      --merge \
      --merged_out "${sample}_merged.fastq" \
      -o "${sample}_unmerged_R1.fastq" \
      -O "${sample}_unmerged_R2.fastq" \
      --unpaired1 "${sample}_singleton_R1.fastq" \
      --unpaired2 "${sample}_singleton_R2.fastq" \
      -h "${sample}_fastp.html" \
      -j "${sample}_fastp.json" \
      --detect_adapter_for_pe \
      --correction \
      --overrepresentation_analysis \
      --length_required 30 \
      --thread ${task.cpus}
    """
}

// ---------------------------------------------------------------------------
// STAGE 2: Dual-Mode Genome Mapping (bwa-meme)
// ---------------------------------------------------------------------------
process BWA_MEME_MAP {
    tag "${sample}"
    publishDir "${params.outdir}/02_mapping_results", mode: 'copy'

    input:
    tuple val(sample), path(merged_fastq), path(unmerged_r1), path(unmerged_r2)
    path ref
    path ref_idx

    output:
    tuple val(sample), path("${sample}_final.sorted.bam"), emit: bam

    script:
    """
    export LC_ALL=C

    # Map Merged Single-End Data
    bwa-meme mem -t ${task.cpus} -R "@RG\\tID:${sample}_merged\\tSM:${sample}\\tPL:ILLUMINA\\tLB:${sample}_lib1" \
    "${ref}" "${merged_fastq}" | \
    samtools sort -@ ${task.cpus} -o "${sample}_merged.sorted.bam" -

    # Map Unmerged Paired-End Data
    bwa-meme mem -t ${task.cpus} -R "@RG\\tID:${sample}_unmerged\\tSM:${sample}\\tPL:ILLUMINA\\tLB:${sample}_lib1" \
    "${ref}" "${unmerged_r1}" "${unmerged_r2}" | \
    samtools sort -@ ${task.cpus} -o "${sample}_unmerged.sorted.bam" -

    # Merge Alignment Streams
    samtools merge -@ ${task.cpus} "${sample}_final.sorted.bam" \
                        "${sample}_merged.sorted.bam" \
                        "${sample}_unmerged.sorted.bam"

    rm "${sample}_merged.sorted.bam" "${sample}_unmerged.sorted.bam"
    """
}

// ---------------------------------------------------------------------------
// STAGE 3: PCR Deduplication (Picard) — failed control batch excluded here
// ---------------------------------------------------------------------------
process MARK_DUPLICATES {
    tag "${sample}"
    publishDir "${params.outdir}/03_deduplicated_bams", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_dedup.bam"), path("${sample}_dedup.bam.bai"), emit: dedup_bam
    path "${sample}_dup_metrics.txt"

    script:
    """
    picard MarkDuplicates \
        I="${bam}" \
        O="${sample}_dedup.bam" \
        M="${sample}_dup_metrics.txt" \
        REMOVE_DUPLICATES=false

    samtools index "${sample}_dedup.bam"
    """
}

// STAGE 3: Multi-Dimensional QC (Picard / Qualimap)
process QC_METRICS {
    tag "${sample}"

    publishDir "${params.outdir}/04_qc_metrics/flagstat", mode: 'copy', pattern: "*_dedup_flagstat.txt"
    publishDir "${params.outdir}/04_qc_metrics/picard", mode: 'copy', pattern: "*_{insert_metrics.txt,insert_histogram.pdf,alignment_metrics.txt}"
    publishDir "${params.outdir}/04_qc_metrics/qualimap", mode: 'copy', pattern: "*_qualimap"

    input:
    tuple val(sample), path(dedup_bam), path(dedup_bai)

    output:
    path "${sample}_dedup_flagstat.txt"
    path "${sample}_insert_metrics.txt"
    path "${sample}_insert_histogram.pdf"
    path "${sample}_alignment_metrics.txt"
    path "${sample}_qualimap", type: 'dir'

    script:
    """
    samtools flagstat "${dedup_bam}" > "${sample}_dedup_flagstat.txt"

    picard CollectInsertSizeMetrics \
        I="${dedup_bam}" O="${sample}_insert_metrics.txt" \
        H="${sample}_insert_histogram.pdf" M=0.5

    picard CollectAlignmentSummaryMetrics \
        I="${dedup_bam}" O="${sample}_alignment_metrics.txt"

    qualimap bamqc -bam "${dedup_bam}" -outdir "${sample}_qualimap" \
        --java-mem-size=4G -nt ${task.cpus}
    """
}

// ---------------------------------------------------------------------------
// STAGE 4: Joint Cohort Variant Genotyping & Hard-Filtering (GATK4)
// ---------------------------------------------------------------------------
process HAPLOTYPE_CALLER {
    tag "${sample}"
    publishDir "${params.outdir}/05_gvcf_outputs", mode: 'copy'

    input:
    tuple val(sample), path(dedup_bam), path(dedup_bai)
    path ref
    path ref_idx

    output:
    tuple path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi"), emit: gvcf

    script:
    """
    # Create required fasta dictionary index if missing
    if [ ! -f \$(basename ${ref} .fna).dict ] && [ ! -f \$(basename ${ref} .fasta).dict ]; then
        picard CreateSequenceDictionary R=${ref}
    fi
    samtools faidx ${ref} 2>/dev/null || true

    gatk HaplotypeCaller \
      -R "${ref}" \
      -I "${dedup_bam}" \
      -O "${sample}.g.vcf.gz" \
      -ERC GVCF \
      -ploidy 2 \
      --native-pair-hmm-threads ${task.cpus}

    gatk ValidateVariants -R "${ref}" -V "${sample}.g.vcf.gz" -gvcf
    """
}

process GENOMICSDB_IMPORT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gvcfs
    path tbis
    path intervals
    path ref
    path ref_idx

    output:
    path "ctropicalis_genomics_db", emit: genomicsdb

    script:
    def gatk_samples = gvcfs.collect { "-V ${it}" }.join(' ')
    """
    rm -rf "./ctropicalis_genomics_db"
    gatk --java-options "-Xmx12G" GenomicsDBImport ${gatk_samples} \
        --genomicsdb-workspace-path "./ctropicalis_genomics_db" -L "${intervals}"
    """
}

process GENOTYPE_GVCFS {
    publishDir "${params.outdir}/06_filtered_variants", mode: 'copy'

    input:
    path genomicsdb
    path ref
    path ref_idx

    output:
    path "raw_variants.vcf.gz", emit: raw_vcf

    script:
    """
    gatk GenotypeGVCFs \
        -R ${ref} \
        -V gendb://${genomicsdb} \
        -O raw_variants.vcf.gz
    """
}

process HARD_FILTER_VARIANTS {
    publishDir "${params.outdir}/06_filtered_variants", mode: 'copy'

    input:
    path raw_vcf
    path ref
    path ref_idx

    output:
    path "filtered_snps.vcf.gz"
    path "filtered_indels.vcf.gz"

    script:
    """
    # Extract & filter SNPs
    gatk SelectVariants -R ${ref} -V ${raw_vcf} -select-type SNP -O raw_snps.vcf.gz
    gatk VariantFiltration -R ${ref} -V raw_snps.vcf.gz \
        -filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "snp_hard_filter" -O filtered_snps.vcf.gz

    # Extract & filter Indels
    gatk SelectVariants -R ${ref} -V ${raw_vcf} -select-type INDEL -O raw_indels.vcf.gz
    gatk VariantFiltration -R ${ref} -V raw_indels.vcf.gz \
        -filter "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "indel_hard_filter" -O filtered_indels.vcf.gz
    """
}

// ---------------------------------------------------------------------------
// STAGE 5: Control-Free Copy Number Profiling (CNVkit)
// ---------------------------------------------------------------------------
process CNVKIT_PROFILE {
    tag "${sample}"
    publishDir "${params.outdir}/07_cnvkit_output/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(dedup_bam), path(dedup_bai)
    path ref
    path gtf

    output:
    path "${sample}.cns"
