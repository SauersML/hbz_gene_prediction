#!/usr/bin/env python3
import os
import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from alphagenome.models import dna_client

ROOT = Path('/Users/user/hbz_gene_prediction')
HBZ_FASTA = Path('/Users/user/Downloads/hbz/HBZ_region_UP005hap1.fasta')
COPIES = ['copy1', 'copy2', 'copy3']
ONTOLOGIES = [
    ('UBERON:0002107', 'Fetal/Embryonic Liver'),
    ('UBERON:0001987', 'Placenta (Embryonic)'),
]


def read_fasta(path: Path):
    name = None
    seq = []
    for line in path.read_text().splitlines():
        if not line:
            continue
        if line.startswith('>'):
            name = line[1:].strip()
            continue
        seq.append(re.sub(r'[^ACGTNacgtn]', '', line))
    return name, ''.join(seq).upper()


def parse_attrs(s: str):
    d = {}
    for p in s.strip().split(';'):
        if not p:
            continue
        if '=' in p:
            k, v = p.split('=', 1)
            d[k] = v
    return d


def parse_models(copy: str):
    out = {}

    def collect(path, key_parent=True):
        models = defaultdict(list)
        strands = {}
        for line in path.read_text().splitlines():
            if not line.strip() or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            s, e = int(f[3]), int(f[4])
            attrs = parse_attrs(f[8])
            gid = attrs.get('Parent', attrs.get('ID', f'cds_{s}_{e}')) if key_parent else attrs.get('ID', f'cds_{s}_{e}')
            models[gid].append((s, e))
            strands[gid] = f[6]
        return [{'id': k, 'cds': sorted(v), 'strand': strands[k]} for k, v in models.items()]

    out['AUGUSTUS'] = collect(ROOT / f'{copy}.augustus.gff3', key_parent=True)
    out['GlimmerHMM'] = collect(ROOT / f'{copy}.glimmer.gff3', key_parent=True)

    return out


def smooth(y, w=21):
    if y is None or len(y) < w:
        return y
    k = np.ones(w) / w
    return np.convolve(y, k, mode='same')


def chunk_prob_percent(p, chunk=20):
    # p: per-base probability array in [0,1]
    p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
    starts, ends, centers, pct = [], [], [], []
    n = len(p)
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        w = p[s:e]
        # Probability at least one true splice site in the chunk.
        prob_any = 1.0 - float(np.prod(1.0 - w))
        starts.append(s + 1)
        ends.append(e)
        centers.append((s + 1 + e) / 2.0)
        pct.append(100.0 * prob_any)
    return np.array(starts), np.array(ends), np.array(centers), np.array(pct)


def pad_seq(full_seq: str):
    if len(full_seq) <= dna_client.SEQUENCE_LENGTH_100KB:
        target = dna_client.SEQUENCE_LENGTH_100KB
    else:
        target = dna_client.SEQUENCE_LENGTH_500KB
    left = (target - len(full_seq)) // 2
    right = target - len(full_seq) - left
    return ('N' * left) + full_seq + ('N' * right), left


def run_for_ontology(model, full_seq, ont):
    seq_padded, left = pad_seq(full_seq)
    out = model.predict_sequence(
        sequence=seq_padded,
        requested_outputs=[
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.CAGE,
            dna_client.OutputType.SPLICE_SITES,
        ],
        ontology_terms=[ont],
    )

    sig = {}

    def crop(td):
        vals = td.values[left:left + len(full_seq), :]
        return vals, td.metadata

    vals, meta = crop(out.rna_seq)
    if vals.shape[1] > 0:
        plus = vals[:, (meta['strand'] == '+').values].mean(axis=1) if (meta['strand'] == '+').any() else np.zeros(vals.shape[0])
        minus = vals[:, (meta['strand'] == '-').values].mean(axis=1) if (meta['strand'] == '-').any() else np.zeros(vals.shape[0])
        sig['rna_plus'] = plus
        sig['rna_minus'] = minus

    vals, meta = crop(out.cage)
    if vals.shape[1] > 0:
        plus = vals[:, (meta['strand'] == '+').values].mean(axis=1) if (meta['strand'] == '+').any() else np.zeros(vals.shape[0])
        minus = vals[:, (meta['strand'] == '-').values].mean(axis=1) if (meta['strand'] == '-').any() else np.zeros(vals.shape[0])
        sig['cage_plus'] = plus
        sig['cage_minus'] = minus

    vals, meta = crop(out.splice_sites)
    if vals.shape[1] > 0:
        # Tissue-agnostic donor/acceptor probabilities on +/- strands.
        # Use max across donor/acceptor channels as per-base "any splice site" chance.
        sig['splice_sites_plus'] = vals[:, (meta['strand'] == '+').values].max(axis=1) if (meta['strand'] == '+').any() else np.zeros(vals.shape[0])
        sig['splice_sites_minus'] = vals[:, (meta['strand'] == '-').values].max(axis=1) if (meta['strand'] == '-').any() else np.zeros(vals.shape[0])

    return sig


def main():
    api_key = os.getenv('ALPHAGENOME_API_KEY')
    if not api_key:
        raise SystemExit('Set ALPHAGENOME_API_KEY')

    _, full_seq = read_fasta(HBZ_FASTA)

    copy_seq = {}
    copy_start = {}
    for c in COPIES:
        _, s = read_fasta(ROOT / f'{c}.fa')
        copy_seq[c] = s
        pos = full_seq.find(s)
        if pos < 0:
            raise RuntimeError(f'Could not map {c} in full sequence')
        copy_start[c] = pos + 1

    model = dna_client.create(api_key)
    ont_signals = {}
    for ont, label in ONTOLOGIES:
        ont_signals[ont] = run_for_ontology(model, full_seq, ont)

    fig = plt.figure(figsize=(19, 10.2), dpi=220)
    gs = fig.add_gridspec(nrows=4, ncols=3, height_ratios=[1.7, 1, 1, 1], hspace=0.28, wspace=0.22)

    pred_colors = {'AUGUSTUS': '#2563eb', 'GlimmerHMM': '#14b8a6'}
    splice_rows = []

    for j, c in enumerate(COPIES):
        L = len(copy_seq[c])
        s0 = copy_start[c] - 1
        e0 = s0 + L
        gmodels = parse_models(c)

        # Row 1 gene models
        ax = fig.add_subplot(gs[0, j])
        ymap = {'AUGUSTUS': 2.2, 'GlimmerHMM': 1.2}
        for pred, y in ymap.items():
            ax.plot([1, L], [y, y], color='#94a3b8', lw=1.4)
            mods = gmodels.get(pred, [])
            lane_n = max(1, len(mods))
            for idx, m in enumerate(mods):
                yy = y - 0.28 + idx * (0.56 / lane_n)
                cds = m['cds']
                for a, b in zip(cds[:-1], cds[1:]):
                    ax.plot([a[1], b[0]], [yy, yy], color='#64748b', lw=1, linestyle=(0, (2, 2)))
                for s, e in cds:
                    ax.add_patch(Rectangle((s, yy - 0.07), e - s + 1, 0.14, facecolor=pred_colors[pred], edgecolor='black', lw=0.3))
        ax.set_xlim(1, L)
        ax.set_ylim(0.7, 2.7)
        ax.set_yticks([2.2, 1.2])
        ax.set_yticklabels(['AUG', 'Glim'])
        if j == 0:
            ax.set_ylabel('Gene model\npredictor')
        ax.grid(axis='x', color='#e5e7eb', lw=0.5)
        ax.set_title(f'{c} | full {copy_start[c]}-{copy_start[c]+L-1}')

        # Rows 2-3 ontology-specific RNA/CAGE evidence with separate y-axes
        for i, (ont, label) in enumerate(ONTOLOGIES, start=1):
            axo = fig.add_subplot(gs[i, j], sharex=ax)
            axo_cage = axo.twinx()
            sig = ont_signals[ont]
            x = np.arange(1, L + 1)
            has_rna = False
            has_cage = False
            rna_handles, rna_labels = [], []
            cage_handles, cage_labels = [], []
            if 'rna_plus' in sig:
                yp = smooth(sig['rna_plus'][s0:e0])
                h = axo.plot(x, yp, color='#0f766e', lw=0.9, label='RNA +')[0]
                rna_handles.append(h); rna_labels.append('RNA +')
                has_rna = True
            if 'rna_minus' in sig:
                ym = smooth(sig['rna_minus'][s0:e0])
                h = axo.plot(x, ym, color='#14b8a6', lw=0.9, label='RNA -')[0]
                rna_handles.append(h); rna_labels.append('RNA -')
                has_rna = True
            if 'cage_plus' in sig:
                cp = smooth(sig['cage_plus'][s0:e0])
                h = axo_cage.plot(x, cp, color='#b45309', lw=0.9, linestyle='--', label='CAGE +')[0]
                cage_handles.append(h); cage_labels.append('CAGE +')
                has_cage = True
            if 'cage_minus' in sig:
                cm = smooth(sig['cage_minus'][s0:e0])
                h = axo_cage.plot(x, cm, color='#f59e0b', lw=0.9, linestyle='--', label='CAGE -')[0]
                cage_handles.append(h); cage_labels.append('CAGE -')
                has_cage = True

            if has_rna:
                axo.tick_params(axis='y', colors='#0f766e')
            if has_cage:
                axo_cage.tick_params(axis='y', colors='#b45309')

            handles = rna_handles + cage_handles
            labels = rna_labels + cage_labels
            if handles:
                axo.legend(handles, labels, frameon=False, fontsize=7, loc='upper right')
            axo.grid(axis='x', color='#e5e7eb', lw=0.4)
            if j == 0:
                short_label = label.replace(' (Embryonic)','').replace('Fetal/Embryonic ','')
                axo.set_ylabel(short_label + '\nRNA-seq\n(normalized read signal)')
                axo_cage.set_ylabel('CAGE\n(normalized read signal)')
            if i == len(ONTOLOGIES):
                axo.set_xlabel(f'{c} local nt')

        # Row 4: SPLICE_SITES in 20 bp chunks as % chance of >=1 splice site.
        axs = fig.add_subplot(gs[3, j], sharex=ax)
        # use first ontology's splice-site prediction (same across ontology by design)
        base_ont = ONTOLOGIES[0][0]
        s = ont_signals[base_ont]
        has_sp = False
        x = np.arange(1, L + 1)
        if 'splice_sites_plus' in s:
            ysp = smooth(s['splice_sites_plus'][s0:e0])
            bstart, bend, bx, bp_plus = chunk_prob_percent(ysp, chunk=20)
            axs.step(bx, bp_plus, where='mid', color='#1d4ed8', lw=1.0, label='SPLICE + (20bp)')
            has_sp = True
        if 'splice_sites_minus' in s:
            ysm = smooth(s['splice_sites_minus'][s0:e0])
            bstart, bend, bx, bp_minus = chunk_prob_percent(ysm, chunk=20)
            axs.step(bx, bp_minus, where='mid', color='#7c3aed', lw=1.0, label='SPLICE - (20bp)')
            has_sp = True
        if has_sp:
            axs.legend(frameon=False, fontsize=7, loc='upper right')
        axs.grid(axis='x', color='#e5e7eb', lw=0.4)
        axs.set_ylim(0, 100)
        if j == 0:
            axs.set_ylabel('SPLICE_SITES\n% chance in 20 bp bin')
        axs.set_xlabel(f'{c} local nt')

        # Save chunked splice probabilities (plus/minus/any) for this copy.
        if 'splice_sites_plus' in s and 'splice_sites_minus' in s:
            _, _, _, bp_plus = chunk_prob_percent(smooth(s['splice_sites_plus'][s0:e0]), chunk=20)
            bstart, bend, _, bp_minus = chunk_prob_percent(smooth(s['splice_sites_minus'][s0:e0]), chunk=20)
            bp_any = 100.0 * (1.0 - (1.0 - bp_plus / 100.0) * (1.0 - bp_minus / 100.0))
            for k in range(len(bstart)):
                splice_rows.append({
                    'copy': c,
                    'bin_start_local': int(bstart[k]),
                    'bin_end_local': int(bend[k]),
                    'splice_plus_pct': float(bp_plus[k]),
                    'splice_minus_pct': float(bp_minus[k]),
                    'splice_any_pct': float(bp_any[k]),
                })

    fig.suptitle('HBZ copies: gene models with fetal/embryonic AlphaGenome RNA/CAGE evidence', fontsize=14)
    out_png = ROOT / 'hbz_gene_models_with_fetal_alphagenome_tracks.png'
    out_svg = ROOT / 'hbz_gene_models_with_fetal_alphagenome_tracks.svg'
    fig.savefig(out_png)
    fig.savefig(out_svg)

    # Write chunked splice probabilities.
    if splice_rows:
        import pandas as pd
        pd.DataFrame(splice_rows).to_csv(ROOT / 'hbz_splice_sites_20bp_percent.tsv', sep='\t', index=False)

    # mapping table
    with (ROOT / 'hbz_copy_fullregion_mapping.tsv').open('w') as f:
        f.write('copy\tfull_start\tfull_end\tcopy_len\n')
        for c in COPIES:
            s = copy_start[c]
            e = s + len(copy_seq[c]) - 1
            f.write(f'{c}\t{s}\t{e}\t{len(copy_seq[c])}\n')

    print('Wrote', out_png)
    print('Wrote', out_svg)
    if splice_rows:
        print('Wrote', ROOT / 'hbz_splice_sites_20bp_percent.tsv')


if __name__ == '__main__':
    main()
