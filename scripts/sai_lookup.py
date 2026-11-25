#!/usr/bin/env python3
"""Adapt transcript-specific SpliceAI scores to the local Broad SAI-10k module."""

import argparse
import csv
import importlib.util
import json


def read_one(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one row in {path}")
    return rows[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", required=True)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--transcript")
    parser.add_argument("--scores")
    parser.add_argument("--chrom")
    parser.add_argument("--pos", type=int)
    parser.add_argument("--ref")
    parser.add_argument("--alt")
    parser.add_argument("--output")
    parser.add_argument("--json")
    args = parser.parse_args()

    spec = importlib.util.spec_from_file_location("clinv1_sai10k", args.module)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load SAI-10k module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    predict = getattr(module, "sai10k_get_transcript_predictions", None)
    if not callable(predict):
        raise ValueError("SAI-10k module lacks sai10k_get_transcript_predictions")
    if args.check:
        print("SAI-10k prediction function available")
        return
    if any(value is None for value in (args.transcript, args.scores, args.chrom, args.pos,
                                        args.ref, args.alt, args.output, args.json)):
        parser.error("all prediction arguments are required")
    transcript = read_one(args.transcript)
    scores = read_one(args.scores)
    def coord(name):
        value = transcript[name]
        return None if value == "NA" else int(value)
    data = {
        "t_id": transcript["transcript"],
        "STRAND": transcript["strand"],
        "EXON_STARTS": [int(v) for v in transcript["exon_starts"].split(",")],
        "EXON_ENDS": [int(v) for v in transcript["exon_ends"].split(",")],
        "CDS_START": coord("cds_start"),
        "CDS_END": coord("cds_end"),
    }
    for key in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
        data[key] = float(scores[key])
    for key in ("DP_AG", "DP_AL", "DP_DG", "DP_DL"):
        data[key] = int(scores[key])
    result = predict(data, args.pos, chrom=args.chrom, ref=args.ref, alt=args.alt, fasta=None)
    if not isinstance(result, dict) or result.get("transcript_id") != data["t_id"] or not isinstance(result.get("aberrations"), list):
        raise ValueError("Invalid SAI-10k result")
    result.pop("_timing_ms", None)
    with open(args.json, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    columns = ("transcript", "status", "effect", "description", "consequence", "size_bp", "frameshift", "affects_coding")
    rows = []
    for aberration in result["aberrations"]:
        if not isinstance(aberration, dict) or not aberration.get("aberration_type"):
            raise ValueError("Invalid predicted aberration")
        desc = aberration.get("description") or {}
        if not isinstance(desc, dict):
            raise ValueError("Invalid predicted aberration description")
        rows.append({"transcript": data["t_id"], "status": "predicted",
                     "effect": aberration["aberration_type"], "description": desc.get("label") or "",
                     "consequence": desc.get("consequence") or "", "size_bp": desc.get("size_bp"),
                     "frameshift": aberration.get("frameshift"), "affects_coding": aberration.get("affects_coding")})
    if not rows:
        rows.append({"transcript": data["t_id"], "status": "no_effect_predicted", "effect": "none"})
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
