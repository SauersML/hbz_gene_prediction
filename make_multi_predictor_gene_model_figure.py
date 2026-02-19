#!/usr/bin/env python3
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SAMPLES = [
    ("copy1", "copy1.fa", "copy1"),
    ("copy2", "copy2.fa", "copy2"),
    ("copy3", "copy3.fa", "copy3"),
    ("HBZ_hg38", "HBZ_GRCh38_NC000016_152644_154503.fa", "HBZ_GRCh38"),
    ("HBZP1_hg38", "HBZP1_GRCh38_NC000016_163067_165205.fa", "HBZP1_GRCh38"),
]
PREDICTORS = [
    ('AUGUSTUS', '.augustus.gff3'),
    ('GlimmerHMM', '.glimmer.gff3'),
]
ANCHOR_SAMPLE = "HBZ_hg38"
ANCHOR_BP = 100
MAX_SHIFT = 80


def fasta_len(path):
    s = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            continue
        s.append(re.sub(r'[^ACGTNacgtn]', '', line))
    return len(''.join(s))


def fasta_seq(path):
    s = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            continue
        s.append(re.sub(r'[^ACGTNacgtn]', '', line))
    return ''.join(s).upper()


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


def best_anchor_shift(ref_seq, qry_seq, anchor_bp=100, max_shift=80):
    ref = ref_seq[:anchor_bp]
    qry = qry_seq[:anchor_bp]
    best = {'shift': 0, 'matches': -1, 'overlap': -1}
    for shift in range(-max_shift, max_shift + 1):
        matches = 0
        overlap = 0
        for i in range(len(ref)):
            j = i - shift
            if 0 <= j < len(qry):
                overlap += 1
                if ref[i] == qry[j]:
                    matches += 1
        if overlap == 0:
            continue
        if (matches > best['matches']) or (matches == best['matches'] and overlap > best['overlap']):
            best = {'shift': shift, 'matches': matches, 'overlap': overlap}
    ident = best['matches'] / best['overlap'] if best['overlap'] else 0.0
    return best['shift'], best['matches'], best['overlap'], ident


# Collect data
records = []
seq_lengths = {label: fasta_len(fasta) for label, fasta, _ in SAMPLES}
seqs = {label: fasta_seq(fasta) for label, fasta, _ in SAMPLES}

if ANCHOR_SAMPLE not in seqs:
    raise ValueError(f"Missing anchor sample {ANCHOR_SAMPLE}")
anchor_seq = seqs[ANCHOR_SAMPLE]
sample_shift = {}
align_rows = []
for label, _, _ in SAMPLES:
    shift, matches, overlap, ident = best_anchor_shift(
        anchor_seq, seqs[label], anchor_bp=ANCHOR_BP, max_shift=MAX_SHIFT
    )
    sample_shift[label] = shift
    align_rows.append(
        {
            'sample': label,
            'anchor_sample': ANCHOR_SAMPLE,
            'anchor_bp': ANCHOR_BP,
            'best_shift_nt': shift,
            'matches': matches,
            'overlap_bp': overlap,
            'identity': ident,
        }
    )
pd.DataFrame(align_rows).to_csv('hbz_anchor_alignment_100bp.tsv', sep='\t', index=False)

all_models = {}  # (predictor, sample) -> models
for pred, suf in PREDICTORS:
    for label, _, stem in SAMPLES:
        gff = Path(f"{stem}{suf}")
        models = parse_models(gff, pred)
        all_models[(pred, label)] = models
        for m in models:
            cds_nt = sum(e - s + 1 for s, e in m['cds'])
            introns = max(0, len(m['cds']) - 1)
            start = m['cds'][0][0]
            end = m['cds'][-1][1]
            records.append({
                'sample': label,
                'predictor': pred,
                'model_id': m['model_id'],
                'strand': m['strand'],
                'start': start,
                'end': end,
                'cds_exons': len(m['cds']),
                'introns': introns,
                'cds_nt': cds_nt,
                'protein_aa_est': max(0, cds_nt // 3 - 1)
            })

summary = pd.DataFrame(records).sort_values(['predictor', 'sample', 'start'])
summary.to_csv('hbz_multi_predictor_summary_with_hbzp1.tsv', sep='\t', index=False)

# Build single-panel figure with all sample/predictor combinations as lanes
predictor_colors = {
    'AUGUSTUS': '#2563eb',
    'GlimmerHMM': '#dc2626',
}
baseline_color = '#94a3b8'
intron_color = '#475569'

max_len = max(seq_lengths.values())
min_x = 1 + min(sample_shift.values())
max_x = max(seq_lengths[s] + sample_shift[s] for s, _, _ in SAMPLES)
fig_h = max(6.0, 0.95 * len(SAMPLES) + 1.8)
fig, ax = plt.subplots(figsize=(14, fig_h), dpi=220)
group_step = 1.0
pair_offset = 0.08  # tight AUGUSTUS/Glimmer pair
lane_h = 0.13

n_groups = len(SAMPLES)
top_group_y = (n_groups - 1) * group_step

for gidx, (label, _, _) in enumerate(SAMPLES):
    y_center = top_group_y - gidx * group_step
    L = seq_lengths[label]
    xshift = sample_shift[label]

    # one shared baseline per sequence group
    ax.plot([1 + xshift, L + xshift], [y_center, y_center], color=baseline_color, lw=2, solid_capstyle='round')

    for pred, _ in PREDICTORS:
        y = y_center + (pair_offset if pred == 'AUGUSTUS' else -pair_offset)
        models = all_models[(pred, label)]
        color = predictor_colors.get(pred, '#2563eb')

        for m in models:
            cds = m['cds']
            if len(cds) > 1:
                for a, b in zip(cds[:-1], cds[1:]):
                    ax.plot([a[1] + xshift, b[0] + xshift], [y, y], color=intron_color, lw=1.0, linestyle=(0, (2, 2)))
            for s, e in cds:
                ax.add_patch(Rectangle((s + xshift, y - lane_h / 2), e - s + 1, lane_h, facecolor=color, edgecolor='black', lw=0.35))

        aa_txt = ''
        hit = summary[(summary['predictor'] == pred) & (summary['sample'] == label)]
        if not hit.empty:
            aa_txt = f" | {int(hit.iloc[0]['protein_aa_est'])} aa"
        ax.text(max_x + 20, y, f'{label} | {pred}{aa_txt} | shift={xshift:+d}', va='center', ha='left', fontsize=8.5)

    # Exon AA labels from AUGUSTUS only (numeric labels only).
    aug_models = all_models[('AUGUSTUS', label)]
    if aug_models:
        aug_cds = aug_models[0]['cds']
        y_label = y_center + pair_offset + lane_h / 2 + 0.06
        for s, e in aug_cds:
            exon_aa = (e - s + 1) // 3
            x_mid = (s + e) / 2 + xshift
            ax.text(x_mid, y_label, f'{exon_aa}', ha='center', va='bottom', fontsize=7.5, color='#111827')

ax.set_xlim(min_x, max_x + 520)
ax.set_ylim(-0.5, top_group_y + 0.5)
ax.set_yticks([])
ax.set_xlabel(f'Anchor-aligned local nt coordinate (100 bp anchor to {ANCHOR_SAMPLE}; no liftover)')
ax.set_title('HBZ gene models in one panel: copy1, copy2, copy3, HBZ hg38, HBZP1 hg38 (AUGUSTUS + GlimmerHMM)')
ax.grid(axis='x', color='#e2e8f0', linewidth=0.6)

# light separators between sequence groups
for k in range(1, len(SAMPLES)):
    sep_y = top_group_y - (k - 0.5) * group_step
    ax.axhline(sep_y, color='#e5e7eb', lw=0.8)

legend_handles = [
    Rectangle((0, 0), 1, 1, facecolor=predictor_colors['AUGUSTUS'], edgecolor='black', lw=0.4, label='AUGUSTUS CDS'),
    Rectangle((0, 0), 1, 1, facecolor=predictor_colors['GlimmerHMM'], edgecolor='black', lw=0.4, label='GlimmerHMM CDS'),
]
ax.legend(handles=legend_handles, loc='lower right', frameon=False, fontsize=9)

fig.tight_layout()
fig.savefig('hbz_all_predictors_gene_models_single_panel_with_hbzp1.png')
fig.savefig('hbz_all_predictors_gene_models_single_panel_with_hbzp1.svg')

print('Wrote hbz_multi_predictor_summary_with_hbzp1.tsv')
print('Wrote hbz_all_predictors_gene_models_single_panel_with_hbzp1.png and .svg')
print('Wrote hbz_anchor_alignment_100bp.tsv')
print('\nModel counts by predictor/sample:')
print(summary.groupby(['predictor', 'sample']).size().rename('n_models').to_string())
