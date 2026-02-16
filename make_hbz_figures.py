#!/usr/bin/env python3
import re
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.family'] = 'DejaVu Sans'

# Parse sequence lengths
seq_len = {}
for fa in sorted(Path('.').glob('copy*.fa')):
    name = fa.stem
    s=[]
    for line in fa.read_text().splitlines():
        if line.startswith('>'):
            continue
        s.append(re.sub(r'[^ACGTNacgtn]','',line))
    seq_len[name]=len(''.join(s))

# Parse AUGUSTUS models
rows=[]
for gff in sorted(Path('.').glob('copy*.augustus.gff3')):
    sample = gff.stem.replace('.augustus','')
    gene=None
    cds=[]
    strand='?'
    for line in gff.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f=line.split('\t')
        feat=f[2]
        s,e=int(f[3]),int(f[4])
        strand=f[6]
        if feat=='gene':
            gene=(s,e)
        if feat=='CDS':
            cds.append((s,e))
    cds=sorted(cds)
    introns=[]
    for i in range(len(cds)-1):
        introns.append((cds[i][1]+1,cds[i+1][0]-1))

    cds_nt = sum(e-s+1 for s,e in cds)
    aa_len = cds_nt//3 - 1  # minus terminal stop

    rows.append({
        'sample': sample,
        'seq_len': seq_len.get(sample),
        'strand': strand,
        'gene_start': gene[0] if gene else None,
        'gene_end': gene[1] if gene else None,
        'cds_exons': len(cds),
        'introns': len(introns),
        'cds_nt': cds_nt,
        'protein_aa': aa_len,
        'cds_coords': cds,
        'intron_coords': introns,
    })

df=pd.DataFrame(rows).sort_values('sample')
df[['sample','seq_len','strand','gene_start','gene_end','cds_exons','introns','cds_nt','protein_aa']].to_csv('hbz_gene_prediction_summary.tsv',sep='\t',index=False)

# Figure 1: gene models
fig,ax=plt.subplots(figsize=(11,4.8),dpi=180)
y_positions={s:i for i,s in enumerate(df['sample'])}
for _,r in df.iterrows():
    y=y_positions[r['sample']]
    L=r['seq_len']
    ax.plot([1,L],[y,y],color='#9aa4b2',lw=2,solid_capstyle='round')
    # gene span
    ax.plot([r['gene_start'],r['gene_end']],[y,y],color='#334155',lw=4,alpha=0.7)
    # introns
    for s,e in r['intron_coords']:
        ax.plot([s,e],[y,y],color='#64748b',lw=2,linestyle=(0,(2,2)))
    # exons as boxes
    for s,e in r['cds_coords']:
        ax.add_patch(plt.Rectangle((s,y-0.16),e-s+1,0.32,color='#2563eb',alpha=0.95))
    ax.text(L+30,y,f"{r['protein_aa']} aa",va='center',fontsize=9,color='#111827')

ax.set_yticks(list(y_positions.values()))
ax.set_yticklabels(list(y_positions.keys()))
ax.set_xlabel('Position in input sequence (nt)')
ax.set_title('AUGUSTUS gene models (human parameters) for HBZ candidate sequences')
ax.set_xlim(1,max(df['seq_len'])+250)
ax.set_ylim(-0.7,len(y_positions)-0.3)
ax.grid(axis='x',color='#e5e7eb',lw=0.6)
fig.tight_layout()
fig.savefig('hbz_gene_models.png')

# Figure 2: ORF/protein lengths
fig2,ax2=plt.subplots(figsize=(8,4.6),dpi=180)
colors=['#0ea5e9','#0284c7','#0369a1']
ax2.bar(df['sample'],df['protein_aa'],color=colors[:len(df)],edgecolor='#0f172a',linewidth=0.8)
for i,(_,r) in enumerate(df.iterrows()):
    ax2.text(i,r['protein_aa']+1,f"{r['protein_aa']} aa",ha='center',fontsize=9)
ax2.set_ylabel('Predicted protein length (aa)')
ax2.set_title('Spliced CDS-derived ORF lengths')
ax2.set_ylim(0,max(df['protein_aa'])+20)
ax2.grid(axis='y',color='#e5e7eb',lw=0.6)
fig2.tight_layout()
fig2.savefig('hbz_orf_lengths.png')

print(df[['sample','seq_len','strand','gene_start','gene_end','cds_exons','introns','cds_nt','protein_aa']].to_string(index=False))
print('\nWrote: hbz_gene_prediction_summary.tsv, hbz_gene_models.png, hbz_orf_lengths.png')
