#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for tool in samtools bedtools; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 1; }
done
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'
  printf 'plus\t0\tchr1\t101\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tXS:A:+\n'
  printf 'minus\t16\tchr1\t101\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tXS:A:-\n'
  printf 'mismatch\t0\tchr1\t101\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tXS:A:-\n'
  printf 'untagged\t0\tchr1\t101\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n'
} > "$tmp/input.sam"
samtools view -b "$tmp/input.sam" > "$tmp/input.bam"
printf 'chr1\t1000\n' > "$tmp/sizes"
bash "$root/scripts/stranded_coverage.sh" "$tmp/input.bam" "$tmp/sizes" "$tmp/plus.bedgraph" "$tmp/minus.bedgraph"
printf 'chr1\t100\t110\t1\nchr1\t160\t170\t1\n' > "$tmp/expected_plus"
printf 'chr1\t100\t110\t-1\nchr1\t160\t170\t-1\n' > "$tmp/expected_minus"
diff -u "$tmp/expected_plus" "$tmp/plus.bedgraph"
diff -u "$tmp/expected_minus" "$tmp/minus.bedgraph"
echo "Stranded coverage tests passed."
