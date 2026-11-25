#!/usr/bin/env bash
set -euo pipefail

for tool in samtools sambamba; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 1; }
done
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT

{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'
  printf 'before\t0\tchr1\t20\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n'
  printf 'boundary\t0\tchr1\t91\t60\t20M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n'
  printf 'splice\t0\tchr1\t110\t60\t10M50N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n'
  printf 'inside\t0\tchr1\t150\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n'
  printf 'after\t0\tchr1\t300\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n'
} > "$tmp/input.sam"
samtools view -b "$tmp/input.sam" > "$tmp/input.bam"
samtools index "$tmp/input.bam"
printf 'chr1\t99\t200\n' > "$tmp/region.bed"

sambamba view -h -f bam -L "$tmp/region.bed" "$tmp/input.bam" > "$tmp/old.bam"
samtools view -b "$tmp/input.bam" chr1:100-200 > "$tmp/new.bam"
samtools view "$tmp/old.bam" | LC_ALL=C sort > "$tmp/old.sam"
samtools view "$tmp/new.bam" | LC_ALL=C sort > "$tmp/new.sam"
diff -u "$tmp/old.sam" "$tmp/new.sam"
cut -f 1 "$tmp/new.sam" | LC_ALL=C sort > "$tmp/actual_names"
printf 'boundary\ninside\nsplice\n' > "$tmp/expected_names"
diff -u "$tmp/expected_names" "$tmp/actual_names"
samtools view -b "$tmp/input.bam" chr1:400-500 > "$tmp/empty.bam"
test "$(samtools view -c "$tmp/empty.bam")" = 0
cp "$tmp/input.bam" "$tmp/unindexed.bam"
if samtools view -b "$tmp/unindexed.bam" chr1:100-200 > "$tmp/unindexed-result.bam" 2> "$tmp/unindexed.err"; then
  echo "Expected indexed region extraction to fail without an index" >&2
  exit 1
fi
echo "Indexed region extraction tests passed."
