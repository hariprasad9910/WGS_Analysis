# End-to-End WGS Variant & Copy-Number Nextflow Pipeline (DSL2)

A reproducible, fault-tolerant Nextflow (DSL2) pipeline for diploid whole-genome sequencing (WGS) cohorts (*Candida tropicalis*): automated reference indexing, paired-end read quality control, BWA-MEM alignment, PCR deduplication, multi-dimensional QC, joint cohort variant calling/hard-filtering, target gene variant extraction, and control-free copy-number profiling.

## Pipeline Architecture & Stages

| Stage | Process Name | Tool(s) | Purpose |
|---|---|---|---|
| **0. Reference Preparation** | `STAGE_0_PREPARE_REFERENCE` | `bwa`, `samtools`, `GATK4` | Generates BWA indices, FASTA index (`.fai`), and sequence dictionary (`.dict`) once and broadcasts them to worker nodes. |
| **1. Read Preprocessing** | `STAGE_1_FASTP` | `fastp` | Performs paired-end adapter autodetection, base correction in overlapping regions, polyG/overrepresentation analysis, and length filtering (>= 30 bp). |
| **2. Genome Alignment** | `STAGE_2_MAPPING` | `bwa`, `samtools` | Maps trimmed paired-end reads via `bwa mem` with sample-specific `@RG` headers and pipes directly into coordinate sorting (`samtools sort`). |
| **3. Deduplication & QC** | `STAGE_3A_DEDUPLICATION`<br>`STAGE_3B_ALIGNMENT_QC` | `Picard`, `Qualimap`, `samtools` | Flags PCR/optical duplicates (`MarkDuplicates`) and profiles insert-size distribution, alignment summary metrics, coverage depth, and GC bias. |
| **4. Joint Genotyping & Hard Filtering** | `STAGE_4A_HAPLOTYPE_CALLER`<br>`STAGE_4B_JOINT_GENOTYPING`<br>`STAGE_4C_INTEGRATED_FILTRATION_CONCAT` | `GATK4`, `bcftools`, `vcftools` | Per-sample diploid GVCF calling (`-ERC GVCF -ploidy 2`), cohort consolidation (`GenomicsDBImport`), joint genotyping (`GenotypeGVCFs`), split SNP/Indel hard-filtering (`QD`, `QUAL`, `SOR`, `FS`, `MQ`), C-based `bcftools concat`, and cohort missingness/heterozygosity stats. |
| **5. Downstream Extraction & CNV Profiling** | `STAGE_5A_TARGET_GENE_EXTRACTION`<br>`STAGE_5B_CNVKIT_PROFILING` | `bcftools`, `CNVkit` | Extracts genotype/depth matrices (`GT`, `DP`, `AD`) across target resistance genes and runs control-free batch copy-number/aneuploidy profiling (`cnvkit.py batch -n`). |

## Repository Contents

```text
.
├── README.md          — Pipeline documentation and usage guide
├── LICENSE            — MIT License
├── install.sh         — Automated installer for Nextflow and isolated Conda environments
├── nextflow.config    — Resource allocation, dynamic OOM retry scaling, and execution profiles (local / slurm)
└── main.nf            — Nextflow DSL2 workflow implementation
