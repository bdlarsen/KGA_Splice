#!/usr/bin/env bash
set -euo pipefail

BASE=/drift/RNA/TotalRNA_Splice
DEFAULT_OUTPUT_BASE=/drift/RNA/drift
IMAGE=kga-splice:clinv1-minimal-candidate
EXPECTED_IMAGE=sha256:7d70cc91864f88abd3a2fce0aa21e6be183b07cdd7836c61c7a4cfdd9850ec85
SCRIPTS="$BASE/Pipeline/ClinV1/scripts"
SAI="$BASE/Tools/SpliceAI-lookup/google_cloud_run_services"

usage() {
  cat <<'EOF'
ClinV1 forsknings-/valideringskandidat. Ikke godkendt til klinisk brug.

Interaktivt: bash run_clinv1.sh
Parametre:   bash run_clinv1.sh --variant VAR --gene GEN --tx TRANSKRIPT \
               --samples ID eller ID_NT,ID_PURO --genotype het|hmz|hemi|unknown \
               --bam /absolut/bammappe [--output /absolut/outputmappe]

Variant er CHR-POS-REF-ALT eller en lokal .vcf/.vcf.gz-fil. For et par er
foerste ID uden puromycin og andet ID med puromycin. Enkeltproeve har ukendt
behandling. Identitet, behandling og genotype skal bekraeftes uafhaengigt.
Output er som standard /drift/RNA/drift; --output vaelger en anden eksisterende
outputmappe. Interaktivt bruges standardmappen ved at trykke Enter.
Tilfoej --yes kun efter denne kontrol for at starte uden bekraeftelsesspoergsmaal.
EOF
}

fail() { printf 'FEJL: %s\n' "$*" >&2; exit 1; }
variant= gene= tx= samples= genotype= bam= output_base="$DEFAULT_OUTPUT_BASE" output_set=no confirm=no
if (($# == 0)); then
  [[ -t 0 ]] || fail "Interaktiv indtastning kraever en terminal. Se --help."
  read -r -p "Variant (CHR-POS-REF-ALT eller VCF-sti): " variant
  read -r -p "Gen: " gene
  read -r -p "Versionsspecifikt transkript: " tx
  read -r -p "Proeve-ID (et ID eller ID_NT,ID_PURO): " samples
  read -r -p "Genotype (het/hmz/hemi/unknown): " genotype
  read -r -p "Absolut sti til BAM-mappe: " bam
  read -r -p "Outputmappe [${DEFAULT_OUTPUT_BASE}]: " chosen_output
  if [[ -n "$chosen_output" ]]; then output_base="$chosen_output"; fi
else
  while (($#)); do
    case "$1" in
      --help|-h) usage; exit 0 ;;
      --yes) [[ "$confirm" == no ]] || fail "Gentaget --yes."; confirm=yes; shift ;;
      --variant|--gene|--tx|--samples|--genotype|--bam|--output)
        key="$1"; (($# >= 2)) || fail "Manglende vaerdi for $key."
        [[ -n "$2" && "$2" != --* ]] || fail "Manglende vaerdi for $key."
        case "$key" in
          --variant) [[ -z "$variant" ]] || fail "Gentaget $key."; variant="$2" ;;
          --gene) [[ -z "$gene" ]] || fail "Gentaget $key."; gene="$2" ;;
          --tx) [[ -z "$tx" ]] || fail "Gentaget $key."; tx="$2" ;;
          --samples) [[ -z "$samples" ]] || fail "Gentaget $key."; samples="$2" ;;
          --genotype) [[ -z "$genotype" ]] || fail "Gentaget $key."; genotype="$2" ;;
          --bam) [[ -z "$bam" ]] || fail "Gentaget $key."; bam="$2" ;;
          --output) [[ "$output_set" == no ]] || fail "Gentaget $key."; output_base="$2"; output_set=yes ;;
        esac
        shift 2 ;;
      *) fail "Ukendt argument: $1. Se --help." ;;
    esac
  done
fi

[[ -n "$variant" && -n "$gene" && -n "$tx" && -n "$samples" && -n "$genotype" && -n "$bam" ]] ||
  fail "Alle seks patientoplysninger skal angives. Se --help."
[[ "$gene" =~ ^[A-Za-z0-9._-]+$ && "$tx" =~ ^[A-Za-z0-9._-]+$ ]] || fail "Ugyldigt gen eller transkript."
[[ "$samples" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*(,[A-Za-z0-9][A-Za-z0-9._-]*)?$ ]] || fail "Ugyldigt proeve-ID."
first="${samples%%,*}"
[[ "$samples" != *,* || "$first" != "${samples#*,}" ]] || fail "De to proeve-ID'er skal vaere forskellige."
case "$genotype" in het|hmz|hemi|unknown) ;; *) fail "Genotype skal vaere het, hmz, hemi eller unknown." ;; esac
[[ "$bam" == /* && -d "$bam" ]] || fail "BAM-mappen skal vaere en eksisterende absolut sti."
bam="$(realpath -e -- "$bam")"
[[ "$output_base" == /* && -d "$output_base" ]] || fail "Outputmappen skal vaere en eksisterende absolut sti: $output_base"
output_base="$(realpath -e -- "$output_base")"
[[ -w "$output_base" ]] || fail "Outputmappen er ikke skrivbar: $output_base"

vcf_mount=()
case "$variant" in
  *.vcf|*.vcf.gz)
    [[ "$variant" == /* && -f "$variant" ]] || fail "VCF skal vaere en eksisterende absolut sti."
    variant="$(realpath -e -- "$variant")"
    vcf_mount=(--mount "type=bind,src=$variant,dst=$variant,readonly") ;;
  *)
    [[ "$variant" =~ ^(chr)?([0-9]+|X|Y|M|MT)-[0-9]+-[ACGTacgt]+-[ACGTacgt]+$ ]] ||
      fail "Variant skal vaere CHR-POS-REF-ALT eller en absolut VCF-sti." ;;
esac
for path in "$bam" "$variant" "$BASE" "$output_base"; do
  [[ "$path" != *,* ]] || fail "Stier med komma understottes ikke af Docker --mount: $path"
done

for path in "$SCRIPTS/ClinV1.R" "$SCRIPTS/clinical_core.R" \
  "$BASE/Refs/references/genome.fa" "$BASE/Refs/references/genome.fa.fai" \
  "$BASE/GRCh38_Reference.gtf" "$SAI/sai10k_predictions.py"; do
  [[ -s "$path" ]] || fail "Noedvendig fil mangler eller er tom: $path"
done
for path in "$BASE/VEP_Refs/homo_sapiens_refseq/115_GRCh38" \
  "$BASE/RefSample/TotalRefs"; do
  [[ -d "$path" ]] || fail "Noedvendig mappe mangler: $path"
done
command -v docker >/dev/null || fail "Docker blev ikke fundet."
actual="$(docker image inspect "$IMAGE" --format '{{.Id}}')" || fail "Kan ikke inspicere image $IMAGE."
[[ "$actual" == "$EXPECTED_IMAGE" ]] || fail "Image-ID matcher ikke den testede kandidat: $actual"

echo "Proever: $samples | gen: $gene | transkript: $tx | genotype: $genotype"
echo "Variant: $variant | BAM-mappe: $bam"
if [[ "$samples" == *,* ]]; then
  echo "Foerste ID er nt, andet ID er puro."
else
  echo "Enkeltproeve: behandling er ukendt."
fi
echo "Filer laeses lokalt; Docker-netvaerk er deaktiveret."
echo "Output oprettes i: $output_base"
echo "Kontroller ledig plads paa / og outputfilsystemet:"
df -h / "$output_base"
if [[ "$confirm" != yes ]]; then
  [[ -t 0 ]] || fail "Bekraeftelse kraever en terminal; brug kun --yes efter kontrol."
  read -r -p "Er identitet, behandling og genotype kontrolleret, og skal analysen startes? (ja/nej): " answer
  [[ "$answer" == ja ]] || fail "Analysen blev ikke startet."
fi

run_dir="$(mktemp -d "$output_base/${first}_${gene}_XXXXXXXX")"
mkdir "$run_dir/tmp"
echo "Projektmappe: $run_dir"
{
  printf 'Image: %s\nImage-ID: %s\n' "$IMAGE" "$actual"
  sha256sum "$SCRIPTS/ClinV1.R" "$SCRIPTS/clinical_core.R"
} > "$run_dir/tmp/launcher_provenance.txt"

if docker run --rm --network none \
  --hostname clinv1-local --add-host clinv1-local:127.0.0.1 \
  --mount "type=bind,src=$bam,dst=$bam,readonly" \
  --mount "type=bind,src=$BASE/Refs/references,dst=$BASE/Refs/references,readonly" \
  --mount "type=bind,src=$BASE/GRCh38_Reference.gtf,dst=$BASE/GRCh38_Reference.gtf,readonly" \
  --mount "type=bind,src=$BASE/VEP_Refs,dst=$BASE/VEP_Refs,readonly" \
  --mount "type=bind,src=$BASE/RefSample/TotalRefs,dst=$BASE/RefSample/TotalRefs,readonly" \
  --mount "type=bind,src=$SAI,dst=$SAI,readonly" \
  --mount "type=bind,src=$run_dir,dst=$run_dir" \
  --mount "type=bind,src=$run_dir/tmp,dst=/tmp" \
  --mount "type=bind,src=$SCRIPTS/ClinV1.R,dst=/opt/clinv1/scripts/ClinV1.R,readonly" \
  --mount "type=bind,src=$SCRIPTS/clinical_core.R,dst=/opt/clinv1/scripts/clinical_core.R,readonly" \
  "${vcf_mount[@]}" -e TMPDIR=/tmp "$actual" \
  Rscript /opt/clinv1/scripts/ClinV1.R \
    -variant "$variant" -gene "$gene" -tx "$tx" -project-dir "$run_dir" \
    -bam "$bam" -s "$samples" -gt "$genotype" \
    -g "$BASE/Refs/references/genome.fa" -gtf "$BASE/GRCh38_Reference.gtf" \
    -vep "$BASE/VEP_Refs" -refs "$BASE/RefSample/TotalRefs" \
    -sai-module "$SAI/sai10k_predictions.py"; then
  echo "Koersel afsluttet: $run_dir"
  echo "Gennemgaa STATUS.txt og logs foer fortolkning."
  echo "IGV: $run_dir/IGV/paired_IGV.xml"
  echo "IGV VSwarehouse: $run_dir/IGV/paired_IGV_vswarehouse.xml"
else
  status=$?
  echo "FEJL: Docker/analysen afsluttede med kode $status. Delvis output er bevaret: $run_dir" >&2
  exit "$status"
fi
