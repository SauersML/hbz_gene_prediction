#!/usr/bin/env python3
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

FASTA = Path("HBZP1_GRCh38_NC000016_163067_165205.fa")
AUG_GFF = Path("HBZP1_GRCh38.augustus.gff3")
GLIM_GFF = Path("HBZP1_GRCh38.glimmer.gff3")
SUMMARY_TSV = Path("HBZP1_GRCh38.multi_predictor_summary.tsv")
PLOT_PNG = Path("HBZP1_GRCh38_all_predictors_gene_models.png")
PLOT_SVG = Path("HBZP1_GRCh38_all_predictors_gene_models.svg")


def fasta_len(path: Path) -> int:
    seq = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            continue
        seq.append(re.sub(r"[^ACGTNacgtn]", "", line))
    return len("".join(seq))


def parse_attrs(s: str) -> dict:
    d = {}
    for part in s.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def parse_models(gff_path: Path) -> list:
    groups = defaultdict(list)
    strands = {}

    for line in gff_path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        feat = f[2]
        if feat != "CDS":
            continue
        s, e = int(f[3]), int(f[4])
        strand = f[6]
        attrs = parse_attrs(f[8])
        gid = attrs.get("Parent", attrs.get("ID", f"cds_{s}_{e}"))
        groups[gid].append((s, e))
        strands[gid] = strand

    models = []
    for gid, cds in groups.items():
        cds = sorted(cds)
        models.append({"model_id": gid, "strand": strands.get(gid, "+"), "cds": cds})
    models.sort(key=lambda m: m["cds"][0][0] if m["cds"] else 10**9)
    return models


def run_cmd(cmd: list, out_path: Path) -> None:
    with out_path.open("w") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)


def main() -> None:
    if not FASTA.exists():
        raise FileNotFoundError(f"Missing FASTA: {FASTA}")

    run_cmd(
        ["augustus", "--species=human", "--gff3=on", str(FASTA)],
        AUG_GFF,
    )

    run_cmd(
        [
            "glimmerhmm",
            str(FASTA),
            "/opt/homebrew/Caskroom/miniforge/base/share/glimmerhmm/trained_dir/human",
            "-g",
        ],
        GLIM_GFF,
    )

    seq_length = fasta_len(FASTA)
    predictor_to_models = {
        "AUGUSTUS": parse_models(AUG_GFF),
        "GlimmerHMM": parse_models(GLIM_GFF),
    }

    rows = []
    for predictor, models in predictor_to_models.items():
        for m in models:
            cds_nt = sum(e - s + 1 for s, e in m["cds"])
            rows.append(
                {
                    "sequence": FASTA.stem,
                    "predictor": predictor,
                    "model_id": m["model_id"],
                    "strand": m["strand"],
                    "start": m["cds"][0][0],
                    "end": m["cds"][-1][1],
                    "cds_exons": len(m["cds"]),
                    "introns": max(0, len(m["cds"]) - 1),
                    "cds_nt": cds_nt,
                    "protein_aa_est": max(0, cds_nt // 3 - 1),
                }
            )

    summary = pd.DataFrame(rows).sort_values(["predictor", "start"])
    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 4.8), dpi=220, sharex=True)
    plus_color = "#1d4ed8"
    minus_color = "#b45309"
    baseline_color = "#94a3b8"
    intron_color = "#475569"

    for i, predictor in enumerate(["AUGUSTUS", "GlimmerHMM"]):
        ax = axes[i]
        ax.plot([1, seq_length], [0, 0], color=baseline_color, lw=2)
        models = predictor_to_models[predictor]
        n = max(1, len(models))
        lane_h = 0.7 / n

        for k, m in enumerate(models):
            y = 0.3 - k * lane_h
            color = plus_color if m["strand"] == "+" else minus_color
            cds = m["cds"]

            if len(cds) > 1:
                for a, b in zip(cds[:-1], cds[1:]):
                    ax.plot([a[1], b[0]], [y, y], color=intron_color, lw=1.2, linestyle=(0, (2, 2)))

            for s, e in cds:
                ax.add_patch(
                    Rectangle(
                        (s, y - 0.05),
                        e - s + 1,
                        0.1,
                        facecolor=color,
                        edgecolor="black",
                        lw=0.5,
                    )
                )

        ax.set_yticks([])
        ax.set_ylim(-0.5, 0.45)
        ax.grid(axis="x", color="#e2e8f0", linewidth=0.6)
        ax.set_ylabel(predictor)

    axes[-1].set_xlabel("Position in HBZP1 hg38 locus sequence (nt)")
    axes[-1].set_xlim(1, seq_length)
    fig.suptitle("Canonical HBZP1 (hg38): gene model predictions", fontsize=13)

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=plus_color, edgecolor="black", lw=0.5, label="CDS (+ strand)"),
        Rectangle((0, 0), 1, 1, facecolor=minus_color, edgecolor="black", lw=0.5, label="CDS (- strand)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])
    fig.savefig(PLOT_PNG)
    fig.savefig(PLOT_SVG)

    print(f"Wrote {SUMMARY_TSV}")
    print(f"Wrote {PLOT_PNG}")
    print(f"Wrote {PLOT_SVG}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
