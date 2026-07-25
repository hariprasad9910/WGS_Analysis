# 5-Stage WGS Variant & Copy-Number Nextflow Pipeline

A Nextflow (DSL2) reimplementation of a 5-stage whole-genome sequencing
pipeline: read cleaning, dual-mode mapping, deduplication/QC, joint cohort
variant genotyping, and control-free copy-number profiling.

## Pipeline Stages

| Stage | Tool(s) | Purpose |
|---|---|---|
| 1. Read Cleaning & Fragment Merging | `fastp` | Merges overlapping paired-end fragments into high-fidelity single-end reads to improve mapping accuracy, while tracking overrepresented sequences. Trims low-quality trailing/leading bases; reads shorter than 30 bp are removed. |
| 2. Dual-Mode Genome Mapping | `bwa-meme`, `samtools` | Maps merged and unmerged fractions separately (avoiding distortion of the aligner's insert-size distribution), then merges the streams under a shared `SM` read-group tag per isolate. |
| 3. PCR Deduplication & Multi-Dimensional QC | `Picard`, `Qualimap`, `samtools` | Isolates PCR duplication artifacts from real biology; profiles insert size, alignment summary, and GC-bias. |
| 4. Joint Cohort Variant Genotyping & Hard-Filtering | `GATK4`, `bcftools` | Per-sample GVCF calling, cohort-wide GenomicsDB consolidation, joint genotyping, and hard filtering on `QD`, `QUAL`, `SOR`, `FS`, `MQ`. |
| 5. Control-Free Copy Number Profiling | `CNVkit` | Diagnoses macro-genomic shifts (e.g. trisomies, gene duplications) using an internally corrected synthetic background matrix — no matched-normal panel required. |

## Repository Contents

```
.
├── README.md      — this file
├── install.sh     — installs Nextflow + tool dependencies (conda)
└── main.nf        — the Nextflow DSL2 pipeline
```

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.10)
- conda / mamba
- Tools: `fastp`, `bwa-meme`, `samtools`, `picard`, `qualimap`, `gatk4`,
  `bcftools`, `cnvkit`

Run `install.sh` to set these up automatically in two conda environments
(`env_main` for Stages 1–4, `env_cnvkit` for Stage 5, kept isolated to avoid
dependency conflicts):

```bash
bash install.sh
conda activate env_main
```

## Expected Input Layout

```
raw_fastq/
├── sample01_R1_001.fastq
├── sample01_R2_001.fastq
├── sample02_R1_001.fastq
├── sample02_R2_001.fastq
└── ...

refs/
├── ctropicalis_ref.fna         # bwa-meme indexed reference
├── ctropicalis_ref.fna.*       # bwa-meme index files
└── genomic_intervals.list      # GATK call intervals

refs_gtf/
└── GCF_000006335.3_ASM633v3_genomic.gff
```

## Usage

```bash
nextflow run main.nf \
  --reads       raw_fastq \
  --ref         refs/ctropicalis_ref.fna \
  --gtf         refs_gtf/GCF_000006335.3_ASM633v3_genomic.gff \
  --intervals   refs/genomic_intervals.list \
  --exclude     mk22 \
  --cnvkit_exec $(conda run -n env_cnvkit which cnvkit.py) \
  --outdir      results
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `--reads` | `raw_fastq` | Directory containing `*_R1_001.fastq` / `*_R2_001.fastq` pairs |
| `--ref` | `refs/ctropicalis_ref.fna` | Reference FASTA (bwa-meme indexed) |
| `--gtf` | `refs_gtf/GCF_000006335.3_ASM633v3_genomic.gff` | Annotation file for CNVkit |
| `--intervals` | `refs/genomic_intervals.list` | GATK call intervals |
| `--exclude` | `mk22` | Sample name pattern excluded from Stage 3 onward |
| `--cnvkit_exec` | `cnvkit.py` | Path to the `cnvkit.py` executable |
| `--outdir` | `.` | Output/publish directory |

## Output Structure

```
01_fastp_clean/
02_mapping_results/
03_deduplicated_bams/
04_qc_metrics/
├── flagstat/
├── picard/
└── qualimap/
05_gvcf_outputs/
06_filtered_variants/
└── population_qc/
07_cnvkit_results/
├── sample_outputs/
└── visual_plots/
```

## License

MIT
