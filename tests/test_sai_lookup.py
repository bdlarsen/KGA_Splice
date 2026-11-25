#!/usr/bin/env python3
"""Exercise the adapter against a locally installed Broad module."""

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit("Usage: python test_sai_lookup.py /path/to/sai10k_predictions.py")
module = Path(sys.argv[1]).resolve(strict=True)
adapter = Path(__file__).resolve().parents[1] / "scripts" / "sai_lookup.py"
subprocess.run([sys.executable, str(adapter), "--module", str(module), "--check"], check=True)
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    empty_module = root / "empty.py"
    empty_module.write_text("", encoding="utf-8")
    invalid = subprocess.run([sys.executable, str(adapter), "--module", str(empty_module), "--check"],
                             capture_output=True, text=True)
    assert invalid.returncode != 0
    assert "lacks sai10k_get_transcript_predictions" in invalid.stderr
    transcript = root / "transcript.tsv"
    scores = root / "scores.tsv"
    transcript.write_text(
        "transcript\tstrand\texon_starts\texon_ends\tcds_start\tcds_end\n"
        "NM_test.1\t+\t32889645,32890559,32893214\t32889804,32890664,32893462\t32890598\t32893462\n",
        encoding="utf-8",
    )
    scores.write_text("DS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\n"
                      "0\t0\t0.8\t0.5\t0\t0\t-1\t0\n", encoding="utf-8")
    output = root / "out.tsv"
    full = root / "out.json"
    subprocess.run([sys.executable, str(adapter), "--module", str(module), "--transcript", str(transcript),
                    "--scores", str(scores), "--chrom", "17", "--pos", "32890664", "--ref", "AC",
                    "--alt", "A", "--output", str(output), "--json", str(full)], check=True)
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["effect"] == "partial_exon_deletion"
    assert rows[0]["frameshift"] == "True"
    assert rows[0]["size_bp"] == "1"
    assert rows[0]["transcript"] == "NM_test.1"
    result = json.loads(full.read_text(encoding="utf-8"))
    assert result["aberrations"][0]["stop_codon_introduced"] is None
    scores.write_text("DS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\n"
                      "0\t0\t0\t0\t0\t0\t0\t0\n", encoding="utf-8")
    subprocess.run([sys.executable, str(adapter), "--module", str(module), "--transcript", str(transcript),
                    "--scores", str(scores), "--chrom", "17", "--pos", "32890664", "--ref", "A",
                    "--alt", "G", "--output", str(output), "--json", str(full)], check=True)
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["status"] == "no_effect_predicted"
print("SAI Lookup synthetic SNV and indel tests passed")
