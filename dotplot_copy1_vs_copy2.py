#!/usr/bin/env python3
import re
from pathlib import Path
import matplotlib.pyplot as plt

COMP = str.maketrans('ACGTNacgtn','TGCANtgcan')

def read_fasta(path):
    seq=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            continue
        seq.append(re.sub(r'[^ACGTNacgtn]','',line))
    return ''.join(seq).upper()

def revcomp(s):
    return s.translate(COMP)[::-1]

def augustus_cds_regions(gff_path):
    cds = []
    for line in Path(gff_path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'CDS':
            s, e = int(f[3]), int(f[4])
            cds.append((s, e))
    return sorted(cds)

def kmer_hits(a, b, k=15, max_occ=80):
    # index kmers in b
    idx = {}
    for j in range(len(b)-k+1):
        km = b[j:j+k]
        if 'N' in km:
            continue
        idx.setdefault(km, []).append(j)

    # filter overly repetitive kmers
    idx = {km:pos for km,pos in idx.items() if len(pos) <= max_occ}

    xs, ys = [], []
    for i in range(len(a)-k+1):
        km = a[i:i+k]
        if 'N' in km:
            continue
        js = idx.get(km)
        if not js:
            continue
        for j in js:
            xs.append(i)
            ys.append(j)
    return xs, ys

def main():
    s1 = read_fasta('copy1.fa')
    s2 = read_fasta('copy2.fa')
    s2rc = revcomp(s2)
    cds1 = augustus_cds_regions('copy1.augustus.gff3')
    cds2 = augustus_cds_regions('copy2.augustus.gff3')

    k=15
    max_occ=80
    x_f, y_f = kmer_hits(s1, s2, k=k, max_occ=max_occ)
    x_r, y_r = kmer_hits(s1, s2rc, k=k, max_occ=max_occ)

    fig, axes = plt.subplots(1,2, figsize=(12,5), dpi=220)

    axes[0].scatter(x_f, y_f, s=0.15, alpha=0.7, color='#0f766e', rasterized=True)
    axes[0].set_title(f'Forward dotplot (k={k})')
    axes[0].set_xlabel('copy1 position (nt)')
    axes[0].set_ylabel('copy2 position (nt)')
    axes[0].set_xlim(0, len(s1))
    axes[0].set_ylim(0, len(s2))
    axes[0].grid(color='#e5e7eb', linewidth=0.5)

    axes[1].scatter(x_r, y_r, s=0.15, alpha=0.7, color='#7c3aed', rasterized=True)
    axes[1].set_title(f'Reverse-complement dotplot (k={k})')
    axes[1].set_xlabel('copy1 position (nt)')
    axes[1].set_ylabel('revcomp(copy2) position (nt)')
    axes[1].set_xlim(0, len(s1))
    axes[1].set_ylim(0, len(s2))
    axes[1].grid(color='#e5e7eb', linewidth=0.5)

    # Highlight AUGUSTUS exon regions on relevant axes:
    # x-axis uses copy1 coordinates; y-axis uses copy2 (or revcomp(copy2) for panel 2).
    x_exon_color = '#f59e0b'
    y_exon_color = '#ef4444'
    for s, e in cds1:
        axes[0].axvspan(s, e, color=x_exon_color, alpha=0.10, lw=0)
        axes[1].axvspan(s, e, color=x_exon_color, alpha=0.10, lw=0)
    for s, e in cds2:
        axes[0].axhspan(s, e, color=y_exon_color, alpha=0.10, lw=0)
        # Transform copy2 forward coords to revcomp axis coords.
        rc_s = len(s2) - e + 1
        rc_e = len(s2) - s + 1
        axes[1].axhspan(rc_s, rc_e, color=y_exon_color, alpha=0.10, lw=0)

    fig.suptitle('Dotplot alignment: copy1 vs copy2', fontsize=13)
    fig.text(0.5, 0.015, 'Orange bands: copy1 AUGUSTUS exons (x-axis) | Red bands: copy2 AUGUSTUS exons (y-axis)', ha='center', fontsize=9)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig('copy1_vs_copy2_dotplot_augustus_exons.png')
    fig.savefig('copy1_vs_copy2_dotplot_augustus_exons.svg')

    # quick quantitative summary by fitting diagonal offset in forward space
    offs = [y-x for x,y in zip(x_f,y_f)] if x_f else []
    if offs:
        # robust mode-like summary via binning
        bins = {}
        for d in offs:
            b = int(round(d/5.0)*5)
            bins[b] = bins.get(b,0)+1
        top = sorted(bins.items(), key=lambda t:t[1], reverse=True)[:5]
    else:
        top = []

    with open('copy1_vs_copy2_dotplot_stats.txt','w') as out:
        out.write(f'len(copy1)={len(s1)}\\n')
        out.write(f'len(copy2)={len(s2)}\\n')
        out.write(f'k={k}, max_occ={max_occ}\\n')
        out.write(f'forward_hits={len(x_f)}\\n')
        out.write(f'reverse_hits={len(x_r)}\\n')
        out.write('top_forward_offset_bins(nt):\\n')
        for b,c in top:
            out.write(f'  offset~{b}: {c} hits\\n')

    print(f'Wrote copy1_vs_copy2_dotplot_augustus_exons.png/.svg and stats; forward hits={len(x_f)}, reverse hits={len(x_r)}')

if __name__ == '__main__':
    main()
