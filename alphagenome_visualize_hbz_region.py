#!/usr/bin/env python3
import os
import re
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from alphagenome.models import dna_client

FASTA = Path('/Users/user/Downloads/hbz/HBZ_region_UP005hap1.fasta')
OUTDIR = Path('alphagenome_hbz_visualizations')
OUTDIR.mkdir(exist_ok=True)

COPY_WINDOWS = [
    ('copy1', 1705, 3635),
    ('copy2', 12466, 14642),
    ('copy3', 22484, 24666),
]
ONTOLOGY = 'UBERON:0002106'  # spleen: broad output coverage in this region

REQ = [
    dna_client.OutputType.RNA_SEQ,
    dna_client.OutputType.CAGE,
    dna_client.OutputType.SPLICE_SITES,
    dna_client.OutputType.SPLICE_SITE_USAGE,
    dna_client.OutputType.SPLICE_JUNCTIONS,
    dna_client.OutputType.DNASE,
    dna_client.OutputType.ATAC,
    dna_client.OutputType.CHIP_HISTONE,
    dna_client.OutputType.CHIP_TF,
]


def read_fasta(path: Path) -> str:
    seq = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            continue
        seq.append(re.sub(r'[^ACGTNacgtn]', '', line))
    return ''.join(seq).upper()


def pad_to_supported(seq: str):
    n = len(seq)
    if n <= dna_client.SEQUENCE_LENGTH_16KB:
        L = dna_client.SEQUENCE_LENGTH_16KB
    elif n <= dna_client.SEQUENCE_LENGTH_100KB:
        L = dna_client.SEQUENCE_LENGTH_100KB
    elif n <= dna_client.SEQUENCE_LENGTH_500KB:
        L = dna_client.SEQUENCE_LENGTH_500KB
    else:
        L = dna_client.SEQUENCE_LENGTH_1MB
    left = (L - n) // 2
    right = L - n - left
    return ('N' * left) + seq + ('N' * right), left, L


def crop_track(td, pad_left: int, orig_len: int):
    # Returns x positions in original FASTA coordinates (1-based), cropped values, metadata
    vals = td.values
    meta = td.metadata.copy()
    res = td.resolution

    n = vals.shape[0]
    starts = np.arange(n) * res
    ends = starts + res
    centers = (starts + ends) / 2.0

    # Keep bins whose center lies within original sequence span in padded coordinates
    keep = (centers >= pad_left) & (centers < pad_left + orig_len)
    vals_c = vals[keep, :]
    centers_c = centers[keep] - pad_left + 1  # convert to 1-based original coords
    return centers_c, vals_c, meta, res


def shade_copy_windows(ax):
    for name, s, e in COPY_WINDOWS:
        ax.axvspan(s, e, color='#ef4444', alpha=0.08, lw=0)
        ax.text((s+e)/2, ax.get_ylim()[1] * 0.92, name, ha='center', va='top', fontsize=8, color='#991b1b')


def mean_by_strand(vals, meta, strand):
    if 'strand' not in meta.columns:
        return None
    idx = meta['strand'].astype(str) == strand
    if idx.sum() == 0:
        return None
    return vals[:, idx.values].mean(axis=1)


def save_metadata(out, outdir: Path):
    for name in ['rna_seq', 'cage', 'splice_sites', 'splice_site_usage', 'splice_junctions', 'dnase', 'atac', 'chip_histone', 'chip_tf']:
        if hasattr(out, name):
            try:
                getattr(out, name).metadata.to_csv(outdir / f'{name}_metadata.tsv', sep='\t', index=False)
            except Exception:
                pass


def main():
    api_key = os.getenv('ALPHAGENOME_API_KEY')
    if not api_key:
        raise SystemExit('Missing ALPHAGENOME_API_KEY')
    if not FASTA.exists():
        raise SystemExit(f'Missing FASTA: {FASTA}')

    seq = read_fasta(FASTA)
    orig_len = len(seq)
    seq_pad, pad_left, model_len = pad_to_supported(seq)

    model = dna_client.create(api_key)
    out = model.predict_sequence(sequence=seq_pad, requested_outputs=REQ, ontology_terms=[ONTOLOGY])

    save_metadata(out, OUTDIR)

    meta = {
        'fasta': str(FASTA),
        'orig_len': orig_len,
        'model_len': model_len,
        'pad_left': pad_left,
        'ontology': ONTOLOGY,
    }
    (OUTDIR / 'run_metadata.json').write_text(json.dumps(meta, indent=2))

    # Figure 1: RNA + CAGE
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), dpi=220, sharex=True)

    x, v, m, _ = crop_track(out.rna_seq, pad_left, orig_len)
    y_plus = mean_by_strand(v, m, '+')
    y_minus = mean_by_strand(v, m, '-')
    if y_plus is not None:
        axes[0].plot(x, y_plus, color='#0f766e', lw=0.9, label='RNA +')
    if y_minus is not None:
        axes[0].plot(x, y_minus, color='#14b8a6', lw=0.9, label='RNA -')
    axes[0].set_ylabel('RNA_SEQ')
    axes[0].legend(frameon=False, fontsize=8, loc='upper right')
    axes[0].grid(axis='x', color='#e5e7eb', lw=0.5)

    x, v, m, _ = crop_track(out.cage, pad_left, orig_len)
    y_plus = mean_by_strand(v, m, '+')
    y_minus = mean_by_strand(v, m, '-')
    if y_plus is not None:
        axes[1].plot(x, y_plus, color='#b45309', lw=0.9, label='CAGE +')
    if y_minus is not None:
        axes[1].plot(x, y_minus, color='#f59e0b', lw=0.9, label='CAGE -')
    axes[1].set_ylabel('CAGE')
    axes[1].legend(frameon=False, fontsize=8, loc='upper right')
    axes[1].grid(axis='x', color='#e5e7eb', lw=0.5)

    for ax in axes:
        shade_copy_windows(ax)
    axes[1].set_xlabel('Position in HBZ_region_UP005hap1.fasta (nt)')
    fig.suptitle(f'AlphaGenome RNA/CAGE predictions ({ONTOLOGY})', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUTDIR / 'hbz_rna_cage.png')
    fig.savefig(OUTDIR / 'hbz_rna_cage.svg')

    # Figure 2: Splicing
    fig2, axes2 = plt.subplots(2, 1, figsize=(13, 7), dpi=220, sharex=True)

    x, v, m, _ = crop_track(out.splice_sites, pad_left, orig_len)
    # donor+/acceptor+/donor-/acceptor- expected
    for col, color in [('donor +', '#1d4ed8'), ('acceptor +', '#2563eb'), ('donor -', '#7c3aed'), ('acceptor -', '#a855f7')]:
        mask = ((m['name'].astype(str).str.lower() == col.split()[0]) & (m['strand'].astype(str) == col.split()[1]))
        if mask.sum() > 0:
            axes2[0].plot(x, v[:, mask.values].mean(axis=1), lw=0.8, label=col, color=color)
    axes2[0].set_ylabel('SPLICE_SITES')
    axes2[0].legend(frameon=False, fontsize=8, loc='upper right', ncol=2)
    axes2[0].grid(axis='x', color='#e5e7eb', lw=0.5)

    x, v, m, _ = crop_track(out.splice_site_usage, pad_left, orig_len)
    if v.shape[1] > 0:
        axes2[1].plot(x, v.mean(axis=1), color='#0ea5e9', lw=0.9, label='mean usage')
    axes2[1].set_ylabel('SPLICE_SITE_USAGE')
    axes2[1].legend(frameon=False, fontsize=8, loc='upper right')
    axes2[1].grid(axis='x', color='#e5e7eb', lw=0.5)

    for ax in axes2:
        shade_copy_windows(ax)
    axes2[1].set_xlabel('Position in HBZ_region_UP005hap1.fasta (nt)')
    fig2.suptitle(f'AlphaGenome splicing predictions ({ONTOLOGY})', fontsize=13)
    fig2.tight_layout(rect=[0, 0, 1, 0.97])
    fig2.savefig(OUTDIR / 'hbz_splicing.png')
    fig2.savefig(OUTDIR / 'hbz_splicing.svg')

    # Figure 3: Accessibility + histone/TF (coarse tracks)
    fig3, axes3 = plt.subplots(3, 1, figsize=(13, 8), dpi=220, sharex=True)

    x, v, m, _ = crop_track(out.dnase, pad_left, orig_len)
    if v.shape[1] > 0:
        axes3[0].plot(x, v.mean(axis=1), color='#059669', lw=0.9, label='DNASE')
    x, v, m, _ = crop_track(out.atac, pad_left, orig_len)
    if v.shape[1] > 0:
        axes3[0].plot(x, v.mean(axis=1), color='#10b981', lw=0.9, label='ATAC', alpha=0.8)
    axes3[0].set_ylabel('Accessibility')
    axes3[0].legend(frameon=False, fontsize=8, loc='upper right')
    axes3[0].grid(axis='x', color='#e5e7eb', lw=0.5)

    x, v, m, res = crop_track(out.chip_histone, pad_left, orig_len)
    if v.shape[1] > 0:
        axes3[1].plot(x, v.mean(axis=1), color='#dc2626', lw=1.1, label=f'CHIP_HISTONE mean (res={res})')
    axes3[1].set_ylabel('Histone ChIP')
    axes3[1].legend(frameon=False, fontsize=8, loc='upper right')
    axes3[1].grid(axis='x', color='#e5e7eb', lw=0.5)

    x, v, m, res = crop_track(out.chip_tf, pad_left, orig_len)
    if v.shape[1] > 0:
        axes3[2].plot(x, v.mean(axis=1), color='#7c2d12', lw=1.1, label=f'CHIP_TF mean (res={res})')
    axes3[2].set_ylabel('TF ChIP')
    axes3[2].legend(frameon=False, fontsize=8, loc='upper right')
    axes3[2].grid(axis='x', color='#e5e7eb', lw=0.5)

    for ax in axes3:
        shade_copy_windows(ax)
    axes3[2].set_xlabel('Position in HBZ_region_UP005hap1.fasta (nt)')
    fig3.suptitle(f'AlphaGenome accessibility/regulatory predictions ({ONTOLOGY})', fontsize=13)
    fig3.tight_layout(rect=[0, 0, 1, 0.97])
    fig3.savefig(OUTDIR / 'hbz_accessibility_regulatory.png')
    fig3.savefig(OUTDIR / 'hbz_accessibility_regulatory.svg')

    # Save quick stats
    stats = {
        'rna_tracks': int(out.rna_seq.values.shape[1]),
        'cage_tracks': int(out.cage.values.shape[1]),
        'splice_sites_tracks': int(out.splice_sites.values.shape[1]),
        'splice_site_usage_tracks': int(out.splice_site_usage.values.shape[1]),
        'splice_junctions_tracks': int(out.splice_junctions.values.shape[1]),
        'dnase_tracks': int(out.dnase.values.shape[1]),
        'atac_tracks': int(out.atac.values.shape[1]),
        'chip_histone_tracks': int(out.chip_histone.values.shape[1]),
        'chip_tf_tracks': int(out.chip_tf.values.shape[1]),
    }
    (OUTDIR / 'track_counts.json').write_text(json.dumps(stats, indent=2))

    print('Generated visualizations in', OUTDIR)
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
