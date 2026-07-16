#!/bin/bash

# Create the output directory
mkdir -p fastp

# Loop through all R1 files
for r1 in *_R1_001.fastq; do
    # Derive the corresponding R2 file name
    r2="${r1/_R1_/_R2_}"
    # Extract the base sample name
    sample=$(basename "$r1" _R1_001.fastq)

    echo "Processing sample: ${sample}"

    # Execute fastp with the corrected merging flag:--merge resolves overlapping regions in paired-end WGS data to reduce redundancy and improve alignment
    fastp \
      -i "$r1" \
      -I "$r2" \
      -o "fastp/${sample}_unmerged_R1.fastq" \
      -O "fastp/${sample}_unmerged_R2.fastq" \
      --merge \
      --merged_out "fastp/${sample}_merged.fastq" \
      -h "fastp/${sample}_fastp.html" \
      -j "fastp/${sample}_fastp.json" \
      --detect_adapter_for_pe \
      --correction \
      --length_required 30 \
      --thread 4
done
