#!/usr/bin/env bash
set -euo pipefail

[[ $# == 1 && -d "$1" ]] || { echo "Provide an empty scratch directory" >&2; exit 1; }
scratch="$(cd "$1" && pwd)"
[[ -z "$(ls -A "$scratch")" ]] || { echo "Scratch directory must be empty" >&2; exit 1; }
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
base="$scratch/base"
output="$scratch/drift"
mkdir -p "$scratch/bin" "$scratch/bams" "$base/Pipeline/ClinV1/scripts" \
  "$base/Tools/SpliceAI-lookup/google_cloud_run_services" \
  "$base/Refs/references" "$base/VEP_Refs/homo_sapiens_refseq/115_GRCh38" \
  "$base/RefSample/TotalRefs" "$output"
for file in "$base/Pipeline/ClinV1/scripts/ClinV1.R" "$base/Pipeline/ClinV1/scripts/clinical_core.R" \
  "$base/Tools/SpliceAI-lookup/google_cloud_run_services/sai10k_predictions.py" \
  "$base/Refs/references/genome.fa" "$base/Refs/references/genome.fa.fai" \
  "$base/GRCh38_Reference.gtf"; do
  printf 'fixture\n' > "$file"
done
sed -e "s|^BASE=/drift/RNA/TotalRNA_Splice$|BASE=$base|" \
  -e "s|^DEFAULT_OUTPUT_BASE=/drift/RNA/drift$|DEFAULT_OUTPUT_BASE=$output|" \
  "$root/run_clinv1.sh" > "$scratch/launcher.sh"
cat > "$scratch/bin/docker" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  image)
    printf '%s\n' "${MOCK_IMAGE:-sha256:7d70cc91864f88abd3a2fce0aa21e6be183b07cdd7836c61c7a4cfdd9850ec85}" ;;
  run)
    printf '%s\n' "$@" > "$MOCK_ARGS"
    exit "${MOCK_STATUS:-0}" ;;
  *) exit 2 ;;
esac
MOCK
cat > "$scratch/bin/df" <<'MOCK'
#!/usr/bin/env bash
echo "Mock disk report"
MOCK
chmod +x "$scratch/bin/docker" "$scratch/bin/df"
export PATH="$scratch/bin:$PATH" MOCK_ARGS="$scratch/docker_args"

args=(--variant 22-28725241-AC-A --gene CHEK2 --tx NM_007194.4 \
  --samples 123,456 --genotype het --bam "$scratch/bams" --yes)
bash "$scratch/launcher.sh" "${args[@]}"
grep -Fx -- '--network' "$MOCK_ARGS" >/dev/null
grep -Fx -- 'none' "$MOCK_ARGS" >/dev/null
grep -Fx -- '123,456' "$MOCK_ARGS" >/dev/null
grep -Fx -- '-project-dir' "$MOCK_ARGS" >/dev/null
run_dir="$(find "$output" -mindepth 1 -maxdepth 1 -type d -name '123_CHEK2_*' -print)"
[[ -d "$run_dir/tmp" && -f "$run_dir/tmp/launcher_provenance.txt" && ! -d "$run_dir/output" ]]
grep -Fx -- "type=bind,src=$scratch/bams,dst=$scratch/bams,readonly" "$MOCK_ARGS" >/dev/null
grep -Fx -- "type=bind,src=$run_dir/tmp,dst=/tmp" "$MOCK_ARGS" >/dev/null
grep -Fx -- 'sha256:7d70cc91864f88abd3a2fce0aa21e6be183b07cdd7836c61c7a4cfdd9850ec85' "$MOCK_ARGS" >/dev/null

printf 'fixture\n' > "$scratch/variant.vcf"
mkdir "$scratch/custom"
bash "$scratch/launcher.sh" --variant "$scratch/variant.vcf" --gene CHEK2 --tx NM_1.1 \
  --samples 789 --genotype unknown --bam "$scratch/bams" --output "$scratch/custom" --yes
grep -Fx -- "type=bind,src=$scratch/variant.vcf,dst=$scratch/variant.vcf,readonly" "$MOCK_ARGS" >/dev/null
grep -Fx -- '789' "$MOCK_ARGS" >/dev/null
custom_dir="$(find "$scratch/custom" -mindepth 1 -maxdepth 1 -type d -name '789_CHEK2_*' -print)"
[[ -d "$custom_dir/tmp" && ! -d "$custom_dir/output" ]]
grep -Fx -- "type=bind,src=$custom_dir,dst=$custom_dir" "$MOCK_ARGS" >/dev/null

if bash "$scratch/launcher.sh" "${args[@]}" --output "$scratch/nonexistent"; then
  echo "Missing output directory should fail" >&2; exit 1
fi
if bash "$scratch/launcher.sh" "${args[@]}" --output "$output" --output "$scratch/custom"; then
  echo "Duplicate output argument should fail" >&2; exit 1
fi

if bash "$scratch/launcher.sh" --variant 22-1-A-G --gene CHEK2 --tx NM_1.1 \
  --samples 123,123 --genotype het --bam "$scratch/bams" --yes; then
  echo "Duplicate samples should fail" >&2; exit 1
fi
if MOCK_IMAGE=sha256:wrong bash "$scratch/launcher.sh" "${args[@]}"; then
  echo "Unexpected image should fail" >&2; exit 1
fi
if MOCK_STATUS=7 bash "$scratch/launcher.sh" "${args[@]}"; then
  echo "Docker failure should propagate" >&2; exit 1
else
  status=$?
  [[ "$status" == 7 ]] || { echo "Wrong Docker exit status: $status" >&2; exit 1; }
fi
echo "Launcher tests passed."
