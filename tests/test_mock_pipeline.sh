#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
Rscript -e 'if (!requireNamespace("xml2", quietly=TRUE)) stop("Install xml2 in the R environment first")'
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$tmp/bin" "$tmp/bam" "$tmp/vep" "$tmp/out" "$tmp/het" "$tmp/fail" "$tmp/fail_check" "$tmp/fail_extract"
cat > "$tmp/bin/conda" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
input="" output="" tool=""
while (($#)); do
  case "$1" in
    vep|spliceai|sashimipy) tool="$1"; shift ;;
    -i|-I) input="$2"; shift 2 ;;
    -o|-O) output="$2"; shift 2 ;;
    *) shift ;;
  esac
done
if [[ "${FAIL_VEP:-0}" == 1 && "$tool" == vep ]]; then echo "mock failure" >&2; exit 9; fi
if [[ "$tool" == sashimipy ]]; then
  if [[ -z "$output" ]]; then echo "sashimipy help"; exit 0; fi
  printf '%%PDF-1.4\n' > "$output"
  exit 0
elif [[ "$tool" == vep ]]; then
  info='CSQ=G|missense_variant|GENE|NM_1.1'
  printf '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|SYMBOL|Feature">\n' > "$output"
else
  info='SpliceAI=G|GENE|0.1|0.0|0.0|0.0|1|0|0|0'
  : > "$output"
fi
awk -v info="$info" 'BEGIN{OFS="\t"} /^#/ {print;next} {$8=info;print}' "$input" >> "$output"
MOCK
cat > "$tmp/bin/samtools" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  faidx) printf '>chr17:150-150\nA\n' ;;
  idxstats) printf 'chr17\t1000\t1000000\t0\n*\t0\t0\t0\n' ;;
  view)
    if [[ "$2" == "-H" ]]; then
      printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr17\tLN:1000\n'
    elif [[ "$2" == "-c" && "$3" == "-d" ]]; then
      printf '2\n'
    elif [[ "$2" == "-h" && "$3" == "-d" ]]; then
      printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr17\tLN:1000\n'
    elif [[ "$2" == "-b" && "$3" == "-" ]]; then
      cat
    elif [[ "$2" == "-b" && "$3" == "$MOCK_BAM_DIR/"*.bam && "$4" == "chr17:1-1000" && $# == 4 ]]; then
      if [[ "${FAIL_EXTRACT:-0}" == 1 ]]; then echo "mock index failure" >&2; exit 9; fi
      printf 'mock BAM\n'
    else
      echo "Unexpected samtools view arguments: $*" >&2; exit 2
    fi ;;
  quickcheck)
    if [[ "${FAIL_INPUT_CHECK:-0}" == 1 && "$2" == "$MOCK_BAM_DIR/1_reads.bam" ]]; then
      echo "mock quickcheck failure" >&2; exit 9
    fi ;;
  depth)
    if [[ "${*: -1}" == *"_puro."* ]]; then depth=0; else depth=12; fi
    printf 'chr17\t150\t%s\n' "$depth" ;;
  *) exit 2 ;;
esac
MOCK
cat > "$tmp/bin/bedtools" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  bamtobed) cat > /dev/null ;;
  genomecov) cat > /dev/null; printf 'chr17\t100\t110\t5\n' ;;
  *) exit 2 ;;
esac
MOCK
cat > "$tmp/bin/gatk" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  IndexFeatureFile) : ;;
  ASEReadCounter)
    while (($#)); do
      if [[ "$1" == "-O" ]]; then output="$2"; shift 2; else shift; fi
    done
    printf 'contig,position,variantID,refAllele,altAllele,refCount,altCount,totalCount\nchr17,150,.,A,G,12,8,20\n' > "$output" ;;
  *) exit 2 ;;
esac
MOCK
cat > "$tmp/bin/sambamba" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  view)
    if [[ "$2" == "-c" && "$3" == "$MOCK_OUTPUT_DIR/"*.bam && $# == 3 ]]; then
      printf '10\n'
    else
      echo "Unexpected whole-file or extraction sambamba view: $*" >&2; exit 2
    fi ;;
  index) touch "$2.bai" ;;
  *) exit 2 ;;
esac
MOCK
cat > "$tmp/bin/regtools" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${*: -1}" != *"_puro."* ]]; then
  printf 'chr17\t140\t210\tJ1\t10\t+\t140\t210\t0\t2\t10,10\t0,60\n'
fi
MOCK
chmod +x "$tmp/bin/conda" "$tmp/bin/samtools" "$tmp/bin/sambamba" "$tmp/bin/regtools" "$tmp/bin/bedtools" "$tmp/bin/gatk"
export PATH="$tmp/bin:$PATH" MOCK_BAM_DIR="$tmp/bam" MOCK_OUTPUT_DIR="$tmp/out/GENE_1_nt_2_puro/BAM"
touch "$tmp/bam/1_reads.bam" "$tmp/bam/1_reads.bam.bai" "$tmp/bam/2_reads.bam" "$tmp/bam/2_reads.bam.bai"
printf '>chr17\nA\n' > "$tmp/genome.fa"
touch "$tmp/genome.dict"
printf 'chr17\t1000\t7\t1\t2\n' > "$tmp/genome.fa.fai"
printf 'chr17\ttest\texon\t101\t150\t.\t+\t.\tgene_id "GENE"; transcript_id "NM_1.1";\nchr17\ttest\texon\t201\t250\t.\t+\t.\tgene_id "GENE"; transcript_id "NM_1.1";\n' > "$tmp/gene.gtf"
arguments=(-variant 17-150-A-G -gene GENE -tx NM_1.1 -bam "$tmp/bam" -s 1,2 -g "$tmp/genome.fa" -gtf "$tmp/gene.gtf" -vep "$tmp/vep")
Rscript "$root/scripts/ClinV1.R" "${arguments[@]}" -p "$tmp/out"
project="$tmp/out/GENE_1_nt_2_puro"
grep -q 'SpliceAI: annotation present' "$project/STATUS.txt"
grep -q 'zero depth at any variant REF base: 2_puro' "$project/STATUS.txt"
grep -q 'selected_transcript' "$project/Junctions/junction_counts.tsv"
test -s "$project/IGV/paired_IGV.xml"
test -s "$project/Coverage/1_nt.plus.bedgraph"
grep -q 'genotype_not_provided' "$project/ASE/allele_counts.tsv"
grep -q 'NM_1.1' "$project/Report/index.html"
grep -q 'DS_AG' "$project/Report/index.html"
grep -q 'per_million_mapped' "$project/Junctions/junction_counts.tsv"
test -s "$project/Report/junction_normalized.svg"
export MOCK_OUTPUT_DIR="$tmp/het/GENE_1_nt_2_puro/BAM"
Rscript "$root/scripts/ClinV1.R" "${arguments[@]}" -p "$tmp/het" -gt het
grep -q 'counts_available_exploratory' "$tmp/het/GENE_1_nt_2_puro/ASE/allele_counts.tsv"
export MOCK_OUTPUT_DIR="$tmp/out/GENE_1_nt_2_puro/BAM"
if FAIL_VEP=1 Rscript "$root/scripts/ClinV1.R" "${arguments[@]}" -p "$tmp/fail"; then
  echo "Expected VEP failure" >&2
  exit 1
fi
grep -q INCOMPLETE "$tmp/fail/GENE_1_nt_2_puro/STATUS.txt"
if FAIL_INPUT_CHECK=1 Rscript "$root/scripts/ClinV1.R" "${arguments[@]}" -p "$tmp/fail_check"; then
  echo "Expected input quickcheck failure" >&2; exit 1
fi
grep -q INCOMPLETE "$tmp/fail_check/GENE_1_nt_2_puro/STATUS.txt"
grep -q 'mock quickcheck failure' "$tmp/fail_check/GENE_1_nt_2_puro/Logs/1_nt_input_quickcheck.err"
if FAIL_EXTRACT=1 Rscript "$root/scripts/ClinV1.R" "${arguments[@]}" -p "$tmp/fail_extract"; then
  echo "Expected indexed extraction failure" >&2; exit 1
fi
grep -q INCOMPLETE "$tmp/fail_extract/GENE_1_nt_2_puro/STATUS.txt"
grep -q 'mock index failure' "$tmp/fail_extract/GENE_1_nt_2_puro/Logs/1_nt_extract.err"
echo "Mock pipeline tests passed; actual tools and biological performance remain unvalidated."
