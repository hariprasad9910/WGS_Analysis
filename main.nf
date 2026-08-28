#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
    ========================================================================
    C. T R O P I C A L I S   W G S   D S L 2   P I P E L I N E
    ========================================================================
    Reads Input        : ${params.reads}
    Reference Fasta    : ${params.ref_fasta}
    Target Genes       : ${params.target_genes}
    Output Directory   : ${params.outdir}
    ========================================================================
"""

// Define input file channels explicitly
ref_fasta_ch    = file(params.ref_fasta)
ref_gtf_ch      = file(params.ref_gtf)
ref_gff_ch      = file(params.ref_gff)
intervals_ch    = file(params.intervals)
target_genes_ch = Channel.fromPath(params.target_genes).collect()

// Pair raw input fastq structures dynamically
raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

// Pipeline Entry point Workflow block
workflow {
    STAGE_1_FASTP(raw_reads_ch)
    
    // Filter out sample "mk22" explicitly downstream of cleaning to match analytical limits
    cleaned_reads_ch = STAGE_1_FASTP.out.fastq_streams
        .filter { sample_id, files -> !sample_id.contains('mk22') }
        
    STAGE_2_MAPPING(cleaned_reads_ch, ref_fasta_ch)
    STAGE_3A_DEDUPLICATION(STAGE_2_MAPPING.out.raw_bam)
    STAGE_3B_ALIGNMENT_QC(STAGE_3A_DEDUPLICATION.out.dedup_bam, ref_fasta_ch, ref_gtf_ch)
    STAGE_4A_HAPLOTYPE_CALLER(STAGE_3A_DEDUPLICATION.out.dedup_bam, ref_fasta_ch)
    
    // Gather all generated individual GVCF tracks into an aggregate processing list array
    all_gvcfs_ch = STAGE_4A_HAPLOTYPE_CALLER.out.gvcf.collect()
    all_tbis_ch  = STAGE_4A_HAPLOTYPE_CALLER.out.tbi.collect()
    
    STAGE_4B_JOINT_GENOTYPING(all_gvcfs_ch, all_tbis_ch, ref_fasta_ch, intervals_ch)
    STAGE_4C_INTEGRATED_FILTRATION_CONCAT(STAGE_4B_JOINT_GENOTYPING.out.raw_cohort_vcf, ref_fasta_ch)
    
    // Branch structural paths down simultaneously into point-variant screens vs macro-copy grids
    STAGE_5A_TARGET_GENE_EXTRACTION(STAGE_4C_INTEGRATED_FILTRATION_CONCAT.out.clean_vcf, ref_fasta_ch, target_genes_ch)
    
    all_dedup_bams_ch = STAGE_3A_DEDUPLICATION.out.dedup_bam.map { sample_id, bam, bai -> bam }.collect()
    STAGE_5B_CNVKIT_PROFILING(all_dedup_bams_ch, ref_fasta_ch, ref_gff_ch)
}

/*
 * STAGE 1: Fastq Preprocessing & Read Merging via fastp
 */
process STAGE_1_FASTP {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_clean_R1.fastq"), path("*_clean_R2.fastq"), path("*_merged.fastq"), emit: fastq_streams
    path "*.{html,json}"
    
    script:
    """
    fastp \\
      -i ${reads[0]} -I ${reads[1]} \\
      --merge \\
      --merged_out "${sample_id}_merged.fastq" \\
      -o "${sample_id}_clean_R1.fastq" -O "${sample_id}_clean_R2.fastq" \\
      --unpaired1 "${sample_id}_singleton_R1.fastq" --unpaired2 "${sample_id}_singleton_R2.fastq" \\
      -h "${sample_id}_fastp.html" -j "${sample_id}_fastp.json" \\
      --detect_adapter_for_pe --correction --overrepresentation_analysis \\
      --length_required 30 --thread ${task.cpus}
    """
}

/*
 * STAGE 2: Heterogeneous Dual-Stream Mapping with Inline Sorting Execution
 */
process STAGE_2_MAPPING {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(clean_r1), path(clean_r2), path(merged_fout)
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_final.sorted.bam"), emit: raw_bam
    
    script:
    """
    # Force runtime validation mapping indices are local
    if [ ! -f "${ref_fasta}.bwt" ]; then bwa-meme index ${ref_fasta}; fi
    
    bwa-meme mem -t ${task.cpus} -R "@RG\\tID:${sample_id}_merged\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib1" \\
        ${ref_fasta} ${merged_fout} | samtools sort -@ 2 -o tmp_merged.sorted.bam -
        
    bwa-meme mem -t ${task.cpus} -R "@RG\\tID:${sample_id}_unmerged\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib1" \\
        ${ref_fasta} ${clean_r1} ${clean_r2} | samtools sort -@ 2 -o tmp_unmerged.sorted.bam -
        
    samtools merge -@ ${task.cpus} "${sample_id}_final.sorted.bam" tmp_merged.sorted.bam tmp_unmerged.sorted.bam
    rm -f tmp_merged.sorted.bam tmp_unmerged.sorted.bam
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
    """
    picard -Xmx${task.memory.toGiga()}G MarkDuplicates \\
        I=${raw_bam} O="${sample_id}_dedup.bam" M="${sample_id}_dup_metrics.txt" REMOVE_DUPLICATES=false
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
    """
    mkdir -p flagstat picard qualimap
    samtools flagstat ${dedup_bam} > "flagstat/${sample_id}_dedup_flagstat.txt"
    
    picard CollectInsertSizeMetrics I=${dedup_bam} O="picard/${sample_id}_insert_metrics.txt" H="picard/${sample_id}_insert_histogram.pdf" M=0.5
    picard CollectAlignmentSummaryMetrics I=${dedup_bam} O="picard/${sample_id}_alignment_metrics.txt"
    
    qualimap bamqc -bam ${dedup_bam} -outdir "qualimap/${sample_id}_qualimap" --java-mem-size=${task.memory.toGiga()}G -nt ${task.cpus}
    """
}

/*
 * STAGE 4A: Haploid/Diploid Variant Discovery Window via GATK4 GVCF
 */
process STAGE_4A_HAPLOTYPE_CALLER {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path ref_fasta
    
    output:
    path "${sample_id}.g.vcf.gz", emit: gvcf
    path "${sample_id}.g.vcf.gz.tbi", emit: tbi
    
    script:
    """
    # Generate references inline inside sandbox if missing
    if [ ! -f "${ref_fasta}.fai" ]; then samtools faidx ${ref_fasta}; fi
    if [ ! -f "${ref_fasta.baseName}.dict" ]; then gatk CreateSequenceDictionary -R ${ref_fasta}; fi
    
    gatk --java-options "-Xmx${task.memory.toGiga()}G" HaplotypeCaller \\
      -R ${ref_fasta} -I ${dedup_bam} -O "${sample_id}.g.vcf.gz" -ERC GVCF -ploidy 2 --native-pair-hmm-threads ${task.cpus}
    """
}

/*
 * STAGE 4B: Population Architecture Merging via Workspace GenomicsDB
 */
process STAGE_4B_JOINT_GENOTYPING {
    input:
    path all_gvcfs
    path all_tbis
    path ref_fasta
    path intervals
    
    output:
    path "cohort_raw_variants.vcf.gz", emit: raw_cohort_vcf
    
    script:
    """
    # Dynamically build argument list flags maps from incoming collected paths array
    gvcf_inputs=\$(ls *.g.vcf.gz | awk '{print "-V " \$1}' | tr '\\n' ' ')
    
    rm -rf tmp_genomics_db
    gatk --java-options "-Xmx${task.memory.toGiga()}G" GenomicsDBImport \$gvcf_inputs --genomicsdb-workspace-path tmp_genomics_db -L ${intervals}
    gatk --java-options "-Xmx${task.memory.toGiga()}G" GenotypeGVCFs -R ${ref_fasta} -V gendb://tmp_genomics_db -O cohort_raw_variants.vcf.gz
    """
}

/*
 * STAGE 4C: Advanced Split-Channel SNP/Indel Filtration & High-Speed C Concat
 */
process STAGE_4C_INTEGRATED_FILTRATION_CONCAT {
    input:
    path raw_cohort_vcf
    path ref_fasta
    
    output:
    path "cohort_clean_pass_variants.vcf.gz", emit: clean_vcf
    path "cohort_clean_pass_variants.vcf.gz.tbi", emit: clean_tbi
    path "population_qc_metrics/*"
    
    script:
    """
    mkdir -p population_qc_metrics
    bcftools index -t ${raw_cohort_vcf}
    
    # Track 1: Split slice and filter SNPs
    gatk --java-options "-Xmx${task.memory.toGiga()}G" SelectVariants -R ${ref_fasta} -V ${raw_cohort_vcf} -select-type SNP -O raw_snps.vcf.gz
    gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantFiltration -R ${ref_fasta} -V raw_snps.vcf.gz \\
        --filter-expression "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0" --filter-name "SNP_HARD_FILTER" --missing-values-evaluate-as-failing true -O filt_snps.vcf.gz
    bcftools index -t filt_snps.vcf.gz
    
    # Track 2: Split slice and filter Indels 
    gatk --java-options "-Xmx${task.memory.toGiga()}G" SelectVariants -R ${ref_fasta} -V ${raw_cohort_vcf} -select-type INDEL -O raw_indels.vcf.gz
    gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantFiltration -R ${ref_fasta} -V raw_indels.vcf.gz \\
        --filter-expression "QD < 2.0 || QUAL < 30.0 || FS > 200.0 || SOR > 10.0" --filter-name "INDEL_HARD_FILTER" --missing-values-evaluate-as-failing true -O filt_indels.vcf.gz
    bcftools index -t filt_indels.vcf.gz
    
    # Perform clean C-based Concat compilation
    bcftools concat -a filt_snps.vcf.gz filt_indels.vcf.gz -O z -o cohort_filtered_variants.vcf.gz
    bcftools index -t cohort_filtered_variants.vcf.gz
    
    # Extract pristine PASS components into the shared metrics folder
    gatk --java-options "-Xmx${task.memory.toGiga()}G" SelectVariants -R ${ref_fasta} -V cohort_filtered_variants.vcf.gz --exclude-filtered true -O cohort_clean_pass_variants.vcf.gz
    bcftools index -t cohort_clean_pass_variants.vcf.gz
    
    # Run structural database checks
    bcftools stats cohort_clean_pass_variants.vcf.gz > "population_qc_metrics/bcftools_summary.stats"
    vcftools --gzvcf cohort_clean_pass_variants.vcf.gz --missing-indv --out "population_qc_metrics/missingness"
    vcftools --gzvcf cohort_clean_pass_variants.vcf.gz --het --out "population_qc_metrics/heterozygosity"
    """
}

/*
 * STAGE 5A: Precise Coordinate Slicing and Excel-Protected Output Generation
 */
process STAGE_5A_TARGET_GENE_EXTRACTION {
    input:
    path clean_vcf
    path ref_fasta
