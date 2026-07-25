#!/bin/bash
# install.sh
# Installs Nextflow and all tool dependencies required by main.nf,
# using two isolated conda environments (env_main + env_cnvkit) to
# avoid python/dependency conflicts with CNVkit.
set -e

echo "=========================================================="
echo " Installing pipeline dependencies"
echo "=========================================================="

# ---------------------------------------------------------------
# 1. Nextflow
# ---------------------------------------------------------------
if ! command -v nextflow &> /dev/null; then
    echo "Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin/
else
    echo "Nextflow already installed: $(nextflow -version | head -n 1)"
fi

# ---------------------------------------------------------------
# 2. Core pipeline environment (Stages 1-4)
#    fastp, bwa-meme, samtools, picard, qualimap, gatk4, bcftools
# ---------------------------------------------------------------
echo "Creating env_main (fastp, bwa-meme, samtools, picard, qualimap, gatk4, bcftools)..."
conda create -y -n env_main -c bioconda -c conda-forge \
    fastp \
    bwa-meme \
    samtools \
    picard \
    qualimap \
    gatk4 \
    bcftools

# ---------------------------------------------------------------
# 3. Isolated CNVkit environment (Stage 5)
# ---------------------------------------------------------------
echo "Creating env_cnvkit (cnvkit)..."
conda create -y -n env_cnvkit -c bioconda -c conda-forge cnvkit

echo "=========================================================="
echo " Installation complete."
echo " Activate the core environment before running the pipeline:"
echo "   conda activate env_main"
echo ""
echo " main.nf calls CNVkit via the --cnvkit_exec parameter; point"
echo " it at the env_cnvkit binary, e.g.:"
echo "   nextflow run main.nf --cnvkit_exec \$(conda run -n env_cnvkit which cnvkit.py)"
echo "=========================================================="
