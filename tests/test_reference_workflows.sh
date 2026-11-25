#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for tool in Rscript samtools regtools; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 1; }
done
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT

for name in one two; do
  {
    printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'
    printf '%s_a\t0\tchr1\t141\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tXS:A:+\n' "$name"
    printf '%s_b\t0\tchr1\t141\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tXS:A:+\n' "$name"
    printf '%s_c\t0\tchr1\t141\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tXS:A:+\n' "$name"
  } > "$tmp/$name.sam"
  samtools view -b "$tmp/$name.sam" > "$tmp/$name.bam"
  samtools index "$tmp/$name.bam"
done
printf '%s\n%s\n' "$tmp/one.bam" "$tmp/two.bam" > "$tmp/bams.txt"
Rscript "$root/scripts/reference_pool.R" --list "$tmp/bams.txt" --output-dir "$tmp/pool" \
  --condition minus --target-reads 100 --threads 1
test "$(samtools view -c "$tmp/pool/reference_minus.bam")" = 6
test -s "$tmp/pool/reference_minus.bam.bai"
test "$(($(wc -l < "$tmp/pool/inputs.tsv") - 1))" = 2
Rscript "$root/scripts/reference_pool.R" --list "$tmp/bams.txt" --output-dir "$tmp/subsample" \
  --condition plus --target-reads 2 --threads 1
test -s "$tmp/subsample/reference_plus.bam.bai"
printf 'sample_id\tperson_id\tcondition\tbam\tqc\ns1\tp1\tminus\t%s\tpass\n' "$tmp/one.bam" > "$tmp/manifest.tsv"
Rscript "$root/scripts/reference_bank.R" init --bank "$tmp/bank" --reference-id GRCh38_test
Rscript "$root/scripts/reference_bank.R" add --bank "$tmp/bank" --manifest "$tmp/manifest.tsv"
printf 'sample_id\tperson_id\tcondition\tbam\tqc\ns2\tp2\tminus\t%s\tpass\n' "$tmp/two.bam" > "$tmp/manifest.tsv"
Rscript "$root/scripts/reference_bank.R" add --bank "$tmp/bank" --manifest "$tmp/manifest.tsv"
Rscript "$root/scripts/reference_bank.R" add --bank "$tmp/bank" --manifest "$tmp/manifest.tsv"
regtools junctions extract -s XS -a 8 -m 50 -M 500000 "$tmp/one.bam" > "$tmp/patient.bed"
test -s "$tmp/patient.bed"
printf 's2\n' > "$tmp/controls.txt"
Rscript "$root/scripts/reference_bank.R" compare --bank "$tmp/bank" --patient p1 --condition minus \
  --junctions "$tmp/patient.bed" --region chr1:100-300 --controls "$tmp/controls.txt" --output-dir "$tmp/compare"
grep -q observed "$tmp/compare/comparison.tsv"
printf 's1\n' > "$tmp/controls.txt"
if Rscript "$root/scripts/reference_bank.R" compare --bank "$tmp/bank" --patient p1 --condition minus \
  --junctions "$tmp/patient.bed" --region chr1:100-300 --controls "$tmp/controls.txt" --output-dir "$tmp/self"; then
  echo "Expected self-control rejection" >&2
  exit 1
fi
echo "Reference workflow integration tests passed."
