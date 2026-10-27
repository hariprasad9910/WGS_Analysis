#!/bin/bash
# install.sh
# Installs Nextflow and all tool dependencies required by main.nf,
# using two isolated conda environments (env_main + env_cnvkit) to
# avoid Python/R dependency conflicts with CNVkit.
set -euo pipefail

echo "=========================================================="
echo " Installing pipeline dependencies"
echo "=========================================================="

# ---------------------------------------------------------------
# 1. Nextflow (Rootless HPC-safe installation fallback)
# ---------------------------------------------------------------
if ! command -v nextflow &> /dev/null; then
    echo "Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mkdir -p "$HOME/.local/bin"
    mv nextflow "$HOME/.local/bin/"
    export PATH="$HOME/.local/bin:$PATH"
else
    echo "Nextflow already installed: $(nextflow -version | head -n 1)"
fi

# ---------------------------------------------------------------
# 2. Core pipeline environment (Stages 0-5A)
#    Added vcftools, bedtools, and replaced bwa-meme with standard bwa
# ---------------------------------------------------------------
echo "Creating env_main (fastp, bwa, samtools, picard, qualimap, gatk4, bcftools, vcftools, bedtools)..."
conda create -y -n env_main -c bioconda -c conda-forge \
    fastp \
    bwa \
    samtools \
    picard \
    qualimap \
    gatk4 \
    bcftools \
    vcftools \
    bedtools

# ---------------------------------------------------------------
# 3. Isolated CNVkit environment (Stage 5B)
# ---------------------------------------------------------------
echo "Creating env_cnvkit (cnvkit)..."
conda create -y -n env_cnvkit -c bioconda -c conda-forge cnvkit

echo "=========================================================="
echo " Installation complete."
echo " Run locally with:"
echo "   conda activate env_main"
echo "   nextflow run main.nf -profile local --cnvkit_exec \$(conda run -n env_cnvkit which cnvkit.py)"
echo ""
echo " Run on Slurm HPC with:"
echo "   nextflow run main.nf -profile slurm --cnvkit_exec \$(conda run -n env_cnvkit which cnvkit.py)"
echo "=========================================================="
