#!/usr/bin/env bash
# Beregn v6.4-lignende XS-/alignment-strandfiltreret dækning.
set -euo pipefail

if [[ $# != 4 ]]; then
    echo "Usage: bash stranded_coverage.sh input.bam genome.sizes plus.bedgraph minus.bedgraph" >&2
    exit 2
fi
bam=$1
sizes=$2
plus=$3
minus=$4

samtools view -h -d XS:+ "$bam" |
    samtools view -b - |
    bedtools bamtobed -i stdin -split |
    bedtools genomecov -i stdin -bg -g "$sizes" -strand + > "$plus"

samtools view -h -d XS:- "$bam" |
    samtools view -b - |
    bedtools bamtobed -i stdin -split |
    bedtools genomecov -i stdin -bg -g "$sizes" -strand - |
    awk 'BEGIN { OFS="\t" } { $4 = -$4; print }' > "$minus"
