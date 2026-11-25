#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

echo "=== Versions and dependencies ==="
python -c 'import socket; hostname = socket.gethostname(); print("Local hostname:", hostname, socket.gethostbyname(hostname))'
/opt/conda/bin/python -c 'import importlib.metadata as m; import sys; matches = [d for d in m.distributions() if d.metadata["Name"].lower().replace("_", "-") == "conda-anaconda-telemetry"]; sys.exit("Unexpected conda-anaconda-telemetry installation") if matches else print("Conda telemetry plugin is absent")'
Rscript --version
samtools --version | sed -n '1p'
sambamba --version
bedtools --version
gatk_output="$(gatk --version 2>&1)"
printf '%s\n' "$gatk_output"
if printf '%s\n' "$gatk_output" | grep -Eq '(^|[[:space:]])(ERROR|SEVERE)([[:space:]]|:|$)'; then
  echo "GATK reported an error during version check" >&2
  exit 1
fi
command -v regtools
conda run --no-capture-output -n KGA_Splice sashimipy --help > /dev/null
conda run --no-capture-output -n vep_env vep --help > /dev/null
conda run --no-capture-output -n spliceai_env spliceai -h > /dev/null
conda run --no-capture-output -n spliceai_env python -c \
  'import tensorflow, spliceai; print("TensorFlow", tensorflow.__version__)'

echo "=== Synthetic tests ==="
for test in test_core.R test_report.R test_sai.R test_igv_snapshot.R \
            test_report_pipeline.R test_batch_preflight.R \
            test_reference_pool.R test_reference_bank.R; do
  Rscript "tests/$test"
done
for test in test_mock_pipeline.sh test_region_extraction.sh \
            test_stranded_coverage.sh test_reference_workflows.sh; do
  bash "tests/$test"
done
echo "Minimal image synthetic tests passed."
