#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
    ========================================================================
    C. T R O P I C A L I S   W G S   D S L 2   P I P E L I N E
    ========================================================================
    Reads Input        : ${params.reads}
    Reference Fasta    : ${params.ref_fasta}
    Target Genes       : ${params.target_genes}
    Exclude Pattern    : ${params.exclude}
    Output Directory   : ${params.outdir}
    ========================================================================
"""

// Define input file channels
ref_fasta_ch    = file(params.ref_fasta, checkIfExists: true)
ref_gtf_ch      = file(params.ref_gtf, checkIfExists: true)
ref_gff_ch      = file(params.ref_gff, checkIfExists: true)
intervals_ch    = file(params.intervals, checkIfExists: true)
target_genes_ch = Channel.fromPath(params.target_genes, checkIfExists: true).collect()

// Pair raw input fastq files dynamically
raw_reads_ch    = Channel.fromFilePairs(params.reads, checkIfExists: true)

workflow {
    // Index reference genome ONCE and broadcast indices to all parallel tasks
    STAGE_0_PREPARE_REFERENCE(ref_fasta_ch)

    STAGE_1_FASTP(raw_reads_ch)

    // Dynamically filter out excluded sample pattern (if specified)
    cleaned_reads_ch = STAGE_1_FASTP.out.fastq_streams
        .filter { sample_id, r1, r2 -> !params.exclude || !sample_id.contains(params.exclude) }

    STAGE_2_MAPPING(cleaned_reads_ch, STAGE_0_PREPARE_REFERENCE.out.ref_bundle)
    STAGE_3A_DEDUPLICATION(STAGE_2_MAPPING.out.raw_bam)
    STAGE_3B_ALIGNMENT_QC(STAGE_3A_DEDUPLICATION.out.dedup_bam, ref_fasta_ch, ref_gtf_ch)
    STAGE_4A_HAPLOTYPE_CALLER(STAGE_3A_DEDUPLICATION.out.dedup_bam, STAGE_0_PREPARE_REFERENCE.out.ref_bundle)

    // Gather per-sample GVCFs and indices for cohort-wide joint genotyping
    all_gvcfs_ch = STAGE_4A_HAPLOTYPE_CALLER.out.gvcf.collect()
    all_tbis_ch  = STAGE_4A_HAPLOTYPE_CALLER.out.tbi.collect()

    STAGE_4B_JOINT_GENOTYPING(all_gvcfs_ch, all_tbis_ch, STAGE_0_PREPARE_REFERENCE.out.ref_bundle, intervals_ch)
    STAGE_4C_INTEGRATED_FILTRATION_CONCAT(STAGE_4B_JOINT_GENOTYPING.out.raw_cohort_vcf, STAGE_0_PREPARE_REFERENCE.out.ref_bundle)

    // Downstream branches: Target Gene Variant Extraction & Control-Free CNVkit Profiling
    STAGE_5A_TARGET_GENE_EXTRACTION(STAGE_4C_INTEGRATED_FILTRATION_CONCAT.out.clean_vcf, STAGE_0_PREPARE_REFERENCE.out.ref_bundle, target_genes_ch)

    all_dedup_bams_ch = STAGE_3A_DEDUPLICATION.out.dedup_bam.map { sample_id, bam, bai -> bam }.collect()
    all_dedup_bais_ch = STAGE_3A_DEDUPLICATION.out.dedup_bam.map { sample_id, bam, bai -> bai }.collect()
    STAGE_5B_CNVKIT_PROFILING(all_dedup_bams_ch, all_dedup_bais_ch, ref_fasta_ch, ref_gff_ch)
}

/*
 * STAGE 0: Build BWA, SAMtools, and GATK Reference Indices Once
 */
process STAGE_0_PREPARE_REFERENCE {
    tag "${ref_fasta.name}"

    input:
    path ref_fasta

    output:
    tuple path(ref_fasta), path("${ref_fasta}.fai"), path("${ref_fasta.baseName}.dict"), path("${ref_fasta}.*"), emit: ref_bundle

    script:
    """
    bwa index ${ref_fasta}
    samtools faidx ${ref_fasta}
    gatk CreateSequenceDictionary -R ${ref_fasta} -O "${ref_fasta.baseName}.dict"
    """
}

/*
 * STAGE 1: Paired-End Fastq Preprocessing via fastp
 */
process STAGE_1_FASTP {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_clean_R1.fastq.gz"), path("${sample_id}_clean_R2.fastq.gz"), emit: fastq_streams
    path "*.{html,json}"

    script:
    """
    fastp \\
      -i ${reads[0]} -I ${reads[1]} \\
      -o "${sample_id}_clean_R1.fastq.gz" \\
      -O "${sample_id}_clean_R2.fastq.gz" \\
      -h "${sample_id}_fastp.html" \\
      -j "${sample_id}_fastp.json" \\
      --detect_adapter_for_pe \\
      --correction \\
      --overrepresentation_analysis \\
      --length_required 30 \\
      --thread ${task.cpus}
    """
}

/*
 * STAGE 2: Paired-End BWA-MEM Alignment & Coordinate Sorting
 */
process STAGE_2_MAPPING {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(clean_r1), path(clean_r2)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict), path(bwa_indices)

    output:
    tuple val(sample_id), path("${sample_id}_final.sorted.bam"), emit: raw_bam

    script:
    def bwa_threads = Math.max(1, task.cpus - 2)
    """
    bwa mem \\
      -t ${bwa_threads} \\
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib1" \\
      "${ref_fasta}" \\
      "${clean_r1}" \\
      "${clean_r2}" | \\
    samtools sort -@ 2 -o "${sample_id}_final.sorted.bam" -
    """
}

/*
 * STAGE 3A: PCR Duplicate Tracking Control via Picard
 */
process STAGE_3A_DEDUPLICATION {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(raw_bam)

    output:
    tuple val(sample_id), path("${sample_id}_dedup.bam"), path("${sample_id}_dedup.bam.bai"), emit: dedup_bam
    path "${sample_id}_dup_metrics.txt"

    script:
    def avail_mem = Math.max(1, (task.memory.toGiga() * 0.8).intValue())
    """
    picard -Xmx${avail_mem}G MarkDuplicates \\
        I=${raw_bam} \\
        O="${sample_id}_dedup.bam" \\
        M="${sample_id}_dup_metrics.txt" \\
        REMOVE_DUPLICATES=false
    samtools index "${sample_id}_dedup.bam"
    """
}

/*
 * STAGE 3B: Post-Deduplication Multi-Dimensional Quality Matrix Checking
 */
process STAGE_3B_ALIGNMENT_QC {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path ref_fasta
    path ref_gtf

    output:
    path "flagstat/*"
    path "picard/*"
    path "qualimap/*"

    script:
    def avail_mem = Math.max(1, (task.memory.toGiga() * 0.8).intValue())
    """
    mkdir -p flagstat picard qualimap
    samtools flagstat ${dedup_bam} > "flagstat/${sample_id}_dedup_flagstat.txt"

    picard -Xmx${avail_mem}G CollectInsertSizeMetrics \\
        I=${dedup_bam} \\
        O="picard/${sample_id}_insert_metrics.txt" \\
        H="picard/${sample_id}_insert_histogram.pdf" \\
        M=0.5

    picard -Xmx${avail_mem}G CollectAlignmentSummaryMetrics \\
        R=${ref_fasta} \\
        I=${dedup_bam} \\
        O="picard/${sample_id}_alignment_metrics.txt"

    qualimap bamqc \\
        -bam ${dedup_bam} \\
        -gff ${ref_gtf} \\
        -outdir "qualimap/${sample_id}_qualimap" \\
        --java-mem-size=${avail_mem}G \\
        -nt ${task.cpus}
    """
}

/*
 * STAGE 4A: Diploid Variant Discovery via GATK4 HaplotypeCaller (GVCF Mode)
 */
process STAGE_4A_HAPLOTYPE_CALLER {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict), path(bwa_indices)

    output:
    path "${sample_id}.g.vcf.gz", emit: gvcf
    path "${sample_id}.g.vcf.gz.tbi", emit: tbi

    script:
    def avail_mem = Math.max(1, (task.memory.toGiga() * 0.8).intValue())
    """
    gatk --java-options "-Xmx${avail_mem}G" HaplotypeCaller \\
      -R ${ref_fasta} \\
      -I ${dedup_bam} \\
      -O "${sample_id}.g.vcf.gz" \\
      -ERC GVCF \\
      -ploidy 2 \\
      --native-pair-hmm-threads ${task.cpus}
    """
}

/*
 * STAGE 4B: Cohort Consolidation & Joint Genotyping via GenomicsDBImport
 */
process STAGE_4B_JOINT_GENOTYPING {
    input:
    path all_gvcfs
    path all_tbis
    tuple path(ref_fasta), path(ref_fai), path(ref_dict), path(bwa_indices)
    path intervals

    output:
    path "cohort_raw_variants.vcf.gz", emit: raw_cohort_vcf

    script:
    def avail_mem = Math.max(1, (task.memory.toGiga() * 0.8).intValue())
    def gvcf_args = all_gvcfs.collect { "-V ${it}" }.join(' ')
    """
    rm -rf tmp_genomics_db
    gatk --java-options "-Xmx${avail_mem}G" GenomicsDBImport \\
        ${gvcf_args} \\
        --genomicsdb-workspace-path tmp_genomics_db \\
        -L ${intervals}

    gatk --java-options "-Xmx${avail_mem}G" GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -V gendb://tmp_genomics_db \\
        -O cohort_raw_variants.vcf.gz
    """
}

/*
 * STAGE 4C: Split-Channel SNP/Indel Hard Filtration & C-Based bcftools Concat
 */
process STAGE_4C_INTEGRATED_FILTRATION_CONCAT {
    input:
    path raw_cohort_vcf
    tuple path(ref_fasta), path(ref_fai), path(ref_dict), path(bwa_indices)

    output:
    path "cohort_clean_pass_variants.vcf.gz", emit: clean_vcf
    path "cohort_clean_pass_variants.vcf.gz.tbi", emit: clean_tbi
    path "population_qc_metrics/*"

    script:
    def avail_mem = Math.max(1, (task.memory.toGiga() * 0.8).intValue())
    """
    mkdir -p population_qc_metrics
    bcftools index -t ${raw_cohort_vcf}

    # Track 1: Split and filter SNPs
    gatk --java-options "-Xmx${avail_mem}G" SelectVariants \\
        -R ${ref_fasta} -V ${raw_cohort_vcf} -select-type SNP -O raw_snps.vcf.gz

    gatk --java-options "-Xmx${avail_mem}G" VariantFiltration \\
        -R ${ref_fasta} -V raw_snps.vcf.gz \\
        --filter-expression "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0" \\
        --filter-name "SNP_HARD_FILTER" \\
        --missing-values-evaluate-as-failing false \\
        -O filt_snps.vcf.gz
    bcftools index -t filt_snps.vcf.gz

    # Track 2: Split and filter Indels
    gatk --java-options "-Xmx${avail_mem}G" SelectVariants \\
        -R ${ref_fasta} -V ${raw_cohort_vcf} -select-type INDEL -O raw_indels.vcf.gz

    gatk --java-options "-Xmx${avail_mem}G" VariantFiltration \\
        -R ${ref_fasta} -V raw_indels.vcf.gz \\
        --filter-expression "QD < 2.0 || QUAL < 30.0 || FS > 200.0 || SOR > 10.0" \\
        --filter-name "INDEL_HARD_FILTER" \\
        --missing-values-evaluate-as-failing false \\
        -O filt_indels.vcf.gz
    bcftools index -t filt_indels.vcf.gz

    # Concatenate filtered SNPs and Indels
    bcftools concat -a filt_snps.vcf.gz filt_indels.vcf.gz -O z -o cohort_filtered_variants.vcf.gz
    bcftools index -t cohort_filtered_variants.vcf.gz

    # Retain PASS variants only
    gatk --java-options "-Xmx${avail_mem}G" SelectVariants \\
        -R ${ref_fasta} -V cohort_filtered_variants.vcf.gz --exclude-filtered true -O cohort_clean_pass_variants.vcf.gz
    bcftools index -t cohort_clean_pass_variants.vcf.gz

    # Population QC metrics
    bcftools stats cohort_clean_pass_variants.vcf.gz > "population_qc_metrics/bcftools_summary.stats"
    vcftools --gzvcf cohort_clean_pass_variants.vcf.gz --missing-indv --out "population_qc_metrics/missingness"
    vcftools --gzvcf cohort_clean_pass_variants.vcf.gz --het --out "population_qc_metrics/heterozygosity"
    """
}

/*
 * STAGE 5A: Target Gene Variant Coordinate Extraction & TSV Summary
 */
process STAGE_5A_TARGET_GENE_EXTRACTION {
    input:
    path clean_vcf
    tuple path(ref_fasta), path(ref_fai), path(ref_dict), path(bwa_indices)
    path target_genes

    output:
    path "target_gene_variants_summary.tsv"

    script:
    """
    bcftools index -t ${clean_vcf}
    bcftools query -H -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL[\\t%GT\\t%DP\\t%AD]\\n' ${clean_vcf} > target_gene_variants_summary.tsv
    """
}

/*
 * STAGE 5B: Control-Free Copy-Number Profiling via CNVkit
 */
process STAGE_5B_CNVKIT_PROFILING {
    input:
    path all_bams
    path all_bais
    path ref_fasta
    path ref_gff

    output:
    path "sample_outputs/*"
    path "visual_plots/*"

    script:
    """
    mkdir -p sample_outputs visual_plots
    ${params.cnvkit_exec} batch ${all_bams} \\
        -n \\
        -f ${ref_fasta} \\
        --annotate ${ref_gff} \\
        -p ${task.cpus} \\
        -d sample_outputs/

    for cnr in sample_outputs/*.cnr; do
        base=\$(basename "\$cnr" .cnr)
        cns="sample_outputs/\${base}.cns"
        ${params.cnvkit_exec} scatter "\$cnr" -s "\$cns" -o "visual_plots/\${base}_scatter.pdf"
        ${params.cnvkit_exec} diagram "\$cnr" -s "\$cns" -o "visual_plots/\${base}_diagram.pdf"
    done
    """
}
