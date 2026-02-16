#!/usr/bin/env python3
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

COPIES = ['copy1', 'copy2', 'copy3']
PREDICTORS = [
    ('AUGUSTUS', '.augustus.gff3'),
    ('GlimmerHMM', '.glimmer.gff3'),
]


def fasta_len(path):
    s = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            continue
        s.append(re.sub(r'[^ACGTNacgtn]', '', line))
    return len(''.join(s))


def parse_attrs(s):
    d = {}
    for part in s.strip().split(';'):
        if not part:
            continue
        if '=' in part:
            k, v = part.split('=', 1)
            d[k] = v
    return d


def parse_models(gff_path, predictor):
    # returns list of models, each model = dict(strand, cds=[(s,e),...], model_id)
    groups = defaultdict(list)
    strands = {}

    for line in Path(gff_path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9:
            continue
        feat = f[2]
        s, e = int(f[3]), int(f[4])
        strand = f[6]
        attrs = parse_attrs(f[8])

        if feat != 'CDS':
            continue

        gid = attrs.get('Parent', attrs.get('ID', f'cds_{s}_{e}'))

        groups[gid].append((s, e))
        strands[gid] = strand

    models = []
    for gid, cds in groups.items():
        cds = sorted(cds, key=lambda x: x[0])
        models.append({'model_id': gid, 'strand': strands.get(gid, '+'), 'cds': cds})

    # sort models by first CDS coordinate
    models.sort(key=lambda m: m['cds'][0][0] if m['cds'] else 10**9)
    return models


# Collect data
records = []
seq_lengths = {c: fasta_len(f'{c}.fa') for c in COPIES}

all_models = {}  # (predictor, copy) -> models
for pred, suf in PREDICTORS:
    for c in COPIES:
        gff = Path(f'{c}{suf}')
        models = parse_models(gff, pred)
        all_models[(pred, c)] = models
        for m in models:
            cds_nt = sum(e - s + 1 for s, e in m['cds'])
            introns = max(0, len(m['cds']) - 1)
            start = m['cds'][0][0]
            end = m['cds'][-1][1]
            records.append({
                'copy': c,
                'predictor': pred,
                'model_id': m['model_id'],
                'strand': m['strand'],
                'start': start,
                'end': end,
                'cds_exons': len(m['cds']),
                'introns': introns,
                'cds_nt': cds_nt,
                'protein_aa_est': max(0, cds_nt // 3 - 1),
            })

summary = pd.DataFrame(records).sort_values(['predictor', 'copy', 'start'])
summary.to_csv('hbz_multi_predictor_summary.tsv', sep='\t', index=False)

# Build figure
fig, axes = plt.subplots(
    nrows=len(PREDICTORS), ncols=len(COPIES), figsize=(15, 8), dpi=200, sharex=False, sharey=False
)

plus_color = '#1d4ed8'
minus_color = '#b45309'
baseline_color = '#94a3b8'
intron_color = '#475569'

for i, (pred, _) in enumerate(PREDICTORS):
    for j, c in enumerate(COPIES):
        ax = axes[i, j]
        L = seq_lengths[c]
        models = all_models[(pred, c)]

        # Baseline genome coordinate line
        ax.plot([1, L], [0, 0], color=baseline_color, lw=2)

        n = max(1, len(models))
        lane_h = 0.8 / n

        for k, m in enumerate(models):
            y = 0.4 - k * lane_h
            color = plus_color if m['strand'] == '+' else minus_color
            cds = m['cds']

            # connect exons with intron line
            if len(cds) > 1:
                for a, b in zip(cds[:-1], cds[1:]):
                    ax.plot([a[1], b[0]], [y, y], color=intron_color, lw=1.2, linestyle=(0, (2, 2)))

            for s, e in cds:
                ax.add_patch(Rectangle((s, y - 0.045), e - s + 1, 0.09, facecolor=color, edgecolor='black', lw=0.4))

            # strand arrow indicator near model start
            arrow_x = cds[0][0] if m['strand'] == '+' else cds[-1][1]
            txt = '>' if m['strand'] == '+' else '<'
            ax.text(arrow_x, y + 0.085, txt, fontsize=7, color=color, ha='center', va='bottom')

        ax.set_xlim(1, L)
        ax.set_ylim(-0.6, 0.55)
        ax.set_yticks([])
        ax.grid(axis='x', color='#e2e8f0', linewidth=0.6)

        if i == 0:
            ax.set_title(f'{c} ({L} nt)', fontsize=11)
        if j == 0:
            ax.set_ylabel(pred, fontsize=11)
        if i == len(PREDICTORS) - 1:
            ax.set_xlabel('nt position', fontsize=9)

fig.suptitle('HBZ region: gene model predictions across predictors and copies', fontsize=14, y=0.98)

# Legend
handles = [
    Rectangle((0,0),1,1,facecolor=plus_color,edgecolor='black',lw=0.4,label='CDS (+ strand)'),
    Rectangle((0,0),1,1,facecolor=minus_color,edgecolor='black',lw=0.4,label='CDS (- strand)'),
]
fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False, fontsize=10)
fig.tight_layout(rect=[0, 0.05, 1, 0.95])
fig.savefig('hbz_all_predictors_gene_models.png')
fig.savefig('hbz_all_predictors_gene_models.svg')

print('Wrote hbz_multi_predictor_summary.tsv')
print('Wrote hbz_all_predictors_gene_models.png and .svg')
print('\nModel counts by predictor/copy:')
print(summary.groupby(['predictor','copy']).size().rename('n_models').to_string())
