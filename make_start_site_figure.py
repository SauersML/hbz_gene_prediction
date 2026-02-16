#!/usr/bin/env python3
import re
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def read_fa(p):
    s=[]
    for l in Path(p).read_text().splitlines():
        if l.startswith('>'): continue
        s.append(re.sub(r'[^ACGTNacgtn]','',l))
    return ''.join(s).upper()


def parse_aug(gff):
    cds=[];start=None
    for l in Path(gff).read_text().splitlines():
        if l.startswith('#') or not l.strip(): continue
        f=l.split('\t')
        if len(f)<9: continue
        if f[2]=='CDS': cds.append((int(f[3]),int(f[4])))
        if f[2]=='start_codon': start=int(f[3])
    cds=sorted(cds)
    intr=[]
    for i in range(len(cds)-1):
        intr.append((cds[i][1]+1,cds[i+1][0]-1))
    return cds,intr,start


def kozak(seq,pos1):
    i=pos1-1
    b3 = seq[i-3] if i-3>=0 else 'N'
    p4 = seq[i+3] if i+3<len(seq) else 'N'
    return (1 if b3 in 'AG' else 0)+(1 if p4=='G' else 0)

copies=['copy1','copy2','copy3']
fig,axs=plt.subplots(3,1,figsize=(12,7),dpi=220)

for ax,c in zip(axs,copies):
    seq=read_fa(f'{c}.fa')
    cds,intr,start=parse_aug(f'{c}.augustus.gff3')
    L=len(seq)

    ax.plot([1,L],[0,0],color='#94a3b8',lw=2)
    for s,e in intr:
        ax.plot([s,e],[0,0],color='#64748b',lw=2,linestyle=(0,(2,2)))
    for s,e in cds:
        ax.add_patch(Rectangle((s,-0.1),e-s+1,0.2,color='#2563eb',alpha=0.9))

    # ATGs around AUGUSTUS start +/-120
    lo=max(1,start-120); hi=min(L-2,start+120)
    for p in range(lo,hi+1):
        if seq[p-1:p+2] != 'ATG':
            continue
        in_exon=any(a<=p<=b for a,b in cds)
        in_intron=any(a<=p<=b for a,b in intr)
        k=kozak(seq,p)
        if p==start:
            color='#dc2626'; marker='*'; size=85; z=5
        elif in_exon:
            color='#16a34a'; marker='o'; size=32; z=4
        elif in_intron:
            color='#a855f7'; marker='x'; size=40; z=3
        else:
            color='#6b7280'; marker='.'; size=22; z=2
        ax.scatter([p],[0.18],c=color,marker=marker,s=size,zorder=z)
        ax.text(p,0.26,f'{p}\nK{k}',ha='center',va='bottom',fontsize=7,color=color)

    ax.set_xlim(1,L)
    ax.set_ylim(-0.25,0.5)
    ax.set_yticks([])
    ax.grid(axis='x',color='#e5e7eb',lw=0.6)
    ax.set_title(f'{c}: AUGUSTUS CDS structure and nearby ATGs (K=Kozak core score 0-2)',fontsize=10)

axs[-1].set_xlabel('Genomic position (nt)')

legend_items=[
    ('AUGUSTUS start ATG','#dc2626','*'),
    ('Other exon ATG','#16a34a','o'),
    ('Intron ATG','#a855f7','x'),
    ('Other nearby ATG','#6b7280','.'),
]
for lbl,col,m in legend_items:
    axs[0].scatter([],[],c=col,marker=m,s=45,label=lbl)
axs[0].legend(loc='upper right',fontsize=8,frameon=False)

fig.tight_layout()
fig.savefig('start_site_assessment.png')
fig.savefig('start_site_assessment.svg')
print('Wrote start_site_assessment.png/.svg')
