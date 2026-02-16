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
OUTDIR = Path('alphagenome_hbz')
OUTDIR.mkdir(exist_ok=True)

HBZ_WINDOWS = [
    ('copy1_region', 1705, 3635),
    ('copy2_region', 12466, 14642),
    ('copy3_region', 22484, 24666),
]

SUPPORTED = [
    dna_client.SEQUENCE_LENGTH_16KB,
    dna_client.SEQUENCE_LENGTH_100KB,
    dna_client.SEQUENCE_LENGTH_500KB,
    dna_client.SEQUENCE_LENGTH_1MB,
]

ONTOLOGY_CANDIDATES = [
    ['UBERON:0000178'],  # blood (likely)
    ['UBERON:0002371'],  # bone marrow (likely)
    ['UBERON:0002107'],  # liver
    ['UBERON:0002048'],  # lung (known example)
    ['UBERON:0000955'],  # brain (known example)
]

REQUESTED_OUTPUTS = [
    dna_client.OutputType.RNA_SEQ,
    dna_client.OutputType.CAGE,
    dna_client.OutputType.PROCAP,
    dna_client.OutputType.SPLICE_SITES,
    dna_client.OutputType.SPLICE_SITE_USAGE,
]


def read_fasta(path: Path) -> str:
    seq = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            continue
        seq.append(re.sub(r'[^ACGTNacgtn]', '', line))
    return ''.join(seq).upper()


def choose_length(n: int) -> int:
    for s in SUPPORTED:
        if n <= s:
            return s
    raise ValueError(f'sequence too long ({n}) for max supported {SUPPORTED[-1]}')


def center_pad(seq: str, target: int) -> tuple[str, int]:
    if len(seq) == target:
        return seq, 0
    left = (target - len(seq)) // 2
    right = target - len(seq) - left
    return ('N' * left) + seq + ('N' * right), left


def summarize_track(track, track_name: str, pad_left: int, orig_len: int):
    vals = track.values
    # Crop back to original sequence span to avoid padded N regions
    vals = vals[pad_left:pad_left+orig_len, :]
    per_pos = vals.mean(axis=1)

    global_mean = float(np.mean(per_pos))
    global_max = float(np.max(per_pos))
    p95 = float(np.percentile(per_pos, 95))

    rows = []
    for name, s, e in HBZ_WINDOWS:
        s0 = max(0, s-1)
        e0 = min(orig_len, e)
        if e0 <= s0:
            continue
        w = per_pos[s0:e0]
        w_mean = float(np.mean(w))
        w_max = float(np.max(w))
        enrich = (w_mean / global_mean) if global_mean > 0 else np.nan
        rows.append({
            'track': track_name,
            'window': name,
            'start': s,
            'end': e,
            'window_mean': w_mean,
            'window_max': w_max,
            'global_mean': global_mean,
            'global_p95': p95,
            'global_max': global_max,
            'mean_enrichment_vs_global': enrich,
        })
    return rows, per_pos


def main():
    api_key = os.getenv('ALPHAGENOME_API_KEY')
    if not api_key:
        raise SystemExit('Missing ALPHAGENOME_API_KEY environment variable.')
    if not FASTA.exists():
        raise SystemExit(f'FASTA not found: {FASTA}')

    seq = read_fasta(FASTA)
    orig_len = len(seq)
    model_len = choose_length(orig_len)
    seq_padded, pad_left = center_pad(seq, model_len)

    model = dna_client.create(api_key)

    output = None
    ontology_used = None
    errors = []
    for terms in ONTOLOGY_CANDIDATES:
        try:
            output = model.predict_sequence(
                sequence=seq_padded,
                requested_outputs=REQUESTED_OUTPUTS,
                ontology_terms=terms,
            )
            ontology_used = terms
            break
        except Exception as e:
            errors.append({'ontology_terms': terms, 'error': str(e)})

    if output is None:
        raise RuntimeError('All ontology term attempts failed: ' + json.dumps(errors, indent=2))

    # Save metadata
    meta = {
        'fasta': str(FASTA),
        'orig_len': orig_len,
        'model_len': model_len,
        'pad_left': pad_left,
        'ontology_used': ontology_used,
    }
    Path(OUTDIR / 'run_metadata.json').write_text(json.dumps(meta, indent=2))

    track_map = {}
    for name in ['rna_seq', 'cage', 'procap', 'splice_sites', 'splice_site_usage']:
        if hasattr(output, name):
            track_map[name] = getattr(output, name)

    all_rows = []
    per_pos = {}
    for tname, tdata in track_map.items():
        rows, sig = summarize_track(tdata, tname, pad_left, orig_len)
        all_rows.extend(rows)
        per_pos[tname] = sig

        # write track metadata table
        try:
            md = tdata.metadata
            md.to_csv(OUTDIR / f'{tname}_metadata.tsv', sep='\t', index=False)
        except Exception:
            pass

    df = pd.DataFrame(all_rows)
    df.to_csv(OUTDIR / 'window_signal_summary.tsv', sep='\t', index=False)

    # Simple transcription call based on RNA/CAGE/PROCAP enrichment in any HBZ window
    calls = []
    for window in [w[0] for w in HBZ_WINDOWS]:
        sub = df[df['window'] == window]
        if sub.empty:
            continue
        flags = {}
        for track in ['rna_seq', 'cage', 'procap']:
            ss = sub[sub['track'] == track]
            if ss.empty:
                flags[track] = False
                continue
            flags[track] = bool(float(ss['mean_enrichment_vs_global'].iloc[0]) >= 1.5)
        positive = sum(flags.values())
        if positive >= 2:
            verdict = 'supports_transcribed_state'
        elif positive == 1:
            verdict = 'weak_or_context_dependent_support'
        else:
            verdict = 'no_clear_support_in_requested_tracks'
        calls.append({'window': window, **flags, 'verdict': verdict})

    pd.DataFrame(calls).to_csv(OUTDIR / 'transcription_calls.tsv', sep='\t', index=False)

    # Plot major tracks
    fig, axes = plt.subplots(3, 1, figsize=(13, 7), dpi=180, sharex=True)
    plot_tracks = [('rna_seq', '#0f766e'), ('cage', '#b45309'), ('procap', '#7c3aed')]

    x = np.arange(1, orig_len + 1)
    for ax, (tname, color) in zip(axes, plot_tracks):
        if tname not in per_pos:
            ax.text(0.5, 0.5, f'{tname} not available', transform=ax.transAxes, ha='center')
            continue
        y = per_pos[tname]
        ax.plot(x, y, color=color, lw=0.8)
        ax.set_ylabel(tname)
        for wname, s, e in HBZ_WINDOWS:
            ax.axvspan(s, e, color='#ef4444', alpha=0.08, lw=0)
            ax.text((s+e)/2, ax.get_ylim()[1]*0.9, wname, ha='center', fontsize=8, color='#991b1b')
        ax.grid(axis='x', color='#e5e7eb', lw=0.4)

    axes[-1].set_xlabel('Position in HBZ_region_UP005hap1.fasta')
    fig.suptitle('AlphaGenome predicted signals across HBZ region (highlighted windows)', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUTDIR / 'alphagenome_hbz_transcription_tracks.png')
    fig.savefig(OUTDIR / 'alphagenome_hbz_transcription_tracks.svg')

    print('AlphaGenome run complete.')
    print('Output dir:', OUTDIR)
    print('Ontology used:', ontology_used)
    print('Files: run_metadata.json, window_signal_summary.tsv, transcription_calls.tsv, figures + metadata TSVs')


if __name__ == '__main__':
    main()
