#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import pairwise2  # type: ignore[import-untyped]


HBZP1_FASTA = Path("HBZP1_GRCh38_NC000016_163067_165205.fa")
TARGETS = {
    "HBZP1_hg38": Path("HBZP1_GRCh38_NC000016_163067_165205.fa"),
    "copy1": Path("copy1.fa"),
    "copy2": Path("copy2.fa"),
    "copy3": Path("copy3.fa"),
}

# Canonical HBZP1 CDS (AUGUSTUS on HBZP1 hg38)
HBZP1_CDS_EXONS: List[Tuple[int, int]] = [(95, 150), (1416, 1620), (1962, 2090)]

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def read_fasta(path: Path) -> str:
    seq = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            continue
        seq.append(re.sub(r"[^ACGTNacgtn]", "", line))
    return "".join(seq).upper()


def translate(nt: str) -> str:
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append(CODON_TABLE.get(nt[i : i + 3], "X"))
    prot = "".join(aa)
    return prot[:-1] if prot.endswith("*") else prot


def ref_to_target_map(ref: str, target: str) -> Dict[int, int]:
    aln = pairwise2.align.globalms(
        ref,
        target,
        2.0,   # match
        -1.0,  # mismatch
        -10.0, # gap open
        -0.5,  # gap extend
        one_alignment_only=True,
    )[0]
    ref_aln, tgt_aln = aln.seqA, aln.seqB

    mapping: Dict[int, int] = {}
    ref_pos = 0
    tgt_pos = 0
    for rb, tb in zip(ref_aln, tgt_aln):
        if rb != "-":
            ref_pos += 1
        if tb != "-":
            tgt_pos += 1
        if rb != "-" and tb != "-":
            mapping[ref_pos] = tgt_pos
    return mapping


def main() -> None:
    ref = read_fasta(HBZP1_FASTA)

    rows = []
    fasta_out = []
    for sample, fasta_path in TARGETS.items():
        target = read_fasta(fasta_path)
        mapping = ref_to_target_map(ref, target)

        cds_bases: List[str] = []
        exon_cov = []
        for exon_idx, (start, end) in enumerate(HBZP1_CDS_EXONS, start=1):
            mapped = 0
            exon_len = end - start + 1
            for rp in range(start, end + 1):
                tp = mapping.get(rp)
                if tp is None:
                    continue
                mapped += 1
                cds_bases.append(target[tp - 1])
            exon_cov.append((mapped, exon_len))

        cds_nt = "".join(cds_bases)
        protein = translate(cds_nt)

        rows.append(
            {
                "sample": sample,
                "cds_nt_projected": len(cds_nt),
                "protein_aa": len(protein),
                "exon1_mapped_bp": exon_cov[0][0],
                "exon1_ref_bp": exon_cov[0][1],
                "exon2_mapped_bp": exon_cov[1][0],
                "exon2_ref_bp": exon_cov[1][1],
                "exon3_mapped_bp": exon_cov[2][0],
                "exon3_ref_bp": exon_cov[2][1],
                "protein_sequence": protein,
            }
        )

        fasta_out.append(f">{sample}|HBZP1_annotation_transfer|aa={len(protein)}")
        for i in range(0, len(protein), 80):
            fasta_out.append(protein[i : i + 80])

    tsv_path = Path("hbzp1_annotation_transfer_translation.tsv")
    fa_path = Path("hbzp1_annotation_transfer_predicted_proteins.fa")

    header = [
        "sample",
        "cds_nt_projected",
        "protein_aa",
        "exon1_mapped_bp",
        "exon1_ref_bp",
        "exon2_mapped_bp",
        "exon2_ref_bp",
        "exon3_mapped_bp",
        "exon3_ref_bp",
        "protein_sequence",
    ]
    with tsv_path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for r in rows:
            handle.write("\t".join(str(r[c]) for c in header) + "\n")

    fa_path.write_text("\n".join(fasta_out) + "\n", encoding="utf-8")

    print(f"Wrote {tsv_path}")
    print(f"Wrote {fa_path}")
    for r in rows:
        print(f"{r['sample']}\taa={r['protein_aa']}\t{r['protein_sequence']}")


if __name__ == "__main__":
    main()
