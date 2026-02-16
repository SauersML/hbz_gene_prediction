#!/usr/bin/env python3
import re
from pathlib import Path
import pandas as pd


def read_fa(path):
    seq=[]
    for l in Path(path).read_text().splitlines():
        if l.startswith('>'):
            continue
        seq.append(re.sub(r'[^ACGTNacgtn]','',l))
    return ''.join(seq).upper()


def parse_augustus(gff):
    cds=[]
    start=None
    tx_score=None
    tx_range=None
    for l in Path(gff).read_text().splitlines():
        if not l.strip() or l.startswith('#'):
            continue
        f=l.split('\t')
        if len(f)<9:
            continue
        feat=f[2]
        s,e=int(f[3]),int(f[4])
        if feat=='transcript':
            tx_range=(s,e)
            try:
                tx_score=float(f[5])
            except:
                tx_score=None
        if feat=='start_codon':
            start=s
        if feat=='CDS':
            cds.append((s,e))
    cds=sorted(cds)
    introns=[]
    for i in range(len(cds)-1):
        introns.append((cds[i][1]+1,cds[i+1][0]-1))
    return {'start':start,'cds':cds,'introns':introns,'tx_score':tx_score,'tx_range':tx_range}


def in_any(pos, intervals):
    for a,b in intervals:
        if a<=pos<=b:
            return True
    return False


def kozak_core(seq, pos1):
    # pos1 = 1-based A of ATG
    i=pos1-1
    b3 = seq[i-3] if i-3>=0 else 'N'
    p4 = seq[i+3] if i+3 < len(seq) else 'N'
    score = (1 if b3 in 'AG' else 0) + (1 if p4=='G' else 0)
    return score, b3, p4


def upstream_metrics(seq, start1, window=300):
    if start1<=1:
        return {'up_len':0,'gc':0.0,'cpg_oe':0.0,'tataaa':0}
    up = seq[max(0,start1-1-window):start1-1]
    n=len(up)
    gc = 100.0*(up.count('G')+up.count('C'))/n if n else 0.0
    c=up.count('C'); g=up.count('G'); cg=up.count('CG')
    cpg_oe = (cg*n)/(c*g) if n and c and g else 0.0
    tataaa = len(list(re.finditer('TATAAA', up)))
    return {'up_len':n,'gc':gc,'cpg_oe':cpg_oe,'tataaa':tataaa}


rows=[]
for copy in ['copy1','copy2','copy3']:
    seq=read_fa(f'{copy}.fa')
    aug=parse_augustus(f'{copy}.augustus.gff3')
    sampled=parse_augustus(f'{copy}.augustus.sampled.gff3') if Path(f'{copy}.augustus.sampled.gff3').exists() else {'tx_score':None}

    start=aug['start']
    kscore,b3,p4 = kozak_core(seq,start)
    up=upstream_metrics(seq,start,window=300)

    # collect all ATGs inside predicted CDS exons (possible translation starts in mature mRNA)
    exon_atgs=[]
    for s,e in aug['cds']:
        for p in range(s, e-1):
            if seq[p-1:p+2]=='ATG':
                sc,bb3,pp4=kozak_core(seq,p)
                exon_atgs.append((p,sc,bb3,pp4))

    exon_atgs_sorted=sorted(exon_atgs,key=lambda t:(-t[1],abs(t[0]-start),t[0]))
    best_exon = exon_atgs_sorted[0] if exon_atgs_sorted else None
    exon_atg_list=';'.join([f"{p}(score={sc},-3={bb3},+4={pp4})" for p,sc,bb3,pp4 in exon_atgs_sorted[:10]])

    # nearby ATG regardless of exon (to show intron confounders)
    lo=max(1,start-120); hi=min(len(seq)-2,start+120)
    nearby=[]
    for p in range(lo,hi+1):
        if seq[p-1:p+2]=='ATG':
            sc,bb3,pp4=kozak_core(seq,p)
            loc='CDS_exon' if in_any(p,aug['cds']) else ('intron' if in_any(p,aug['introns']) else 'other')
            nearby.append((p,sc,bb3,pp4,loc))
    nearby_sorted=sorted(nearby,key=lambda t:(-t[1],abs(t[0]-start),t[0]))
    nearby_list=';'.join([f"{p}({loc},score={sc})" for p,sc,bb3,pp4,loc in nearby_sorted])

    rows.append({
        'copy':copy,
        'augustus_start':start,
        'augustus_start_kozak_score_0to2':kscore,
        'augustus_-3':b3,
        'augustus_+4':p4,
        'augustus_sampled_tx_score':sampled.get('tx_score'),
        'n_cds_exons':len(aug['cds']),
        'cds_exons':','.join([f'{a}-{b}' for a,b in aug['cds']]),
        'introns':','.join([f'{a}-{b}' for a,b in aug['introns']]),
        'best_exon_atg': f"{best_exon[0]}(score={best_exon[1]})" if best_exon else 'none',
        'all_exon_atgs': exon_atg_list if exon_atg_list else 'none',
        'nearby_atgs_w_location': nearby_list if nearby_list else 'none',
        'upstream_len':up['up_len'],
        'upstream_gc_pct':round(up['gc'],2),
        'upstream_cpg_oe':round(up['cpg_oe'],3),
        'upstream_tataaa_hits':up['tataaa'],
    })


df=pd.DataFrame(rows)
df.to_csv('start_site_assessment.tsv',sep='\t',index=False)

# concise human summary
with open('start_site_assessment_summary.txt','w') as out:
    for _,r in df.iterrows():
        out.write(f"{r['copy']}: AUGUSTUS start={r['augustus_start']} (Kozak {r['augustus_start_kozak_score_0to2']}/2, -3={r['augustus_-3']}, +4={r['augustus_+4']}), sampled_score={r['augustus_sampled_tx_score']}\\n")
        out.write(f"  CDS exons: {r['cds_exons']}\\n")
        out.write(f"  Introns: {r['introns']}\\n")
        out.write(f"  Exon-contained ATGs: {r['all_exon_atgs']}\\n")
        out.write(f"  Nearby ATGs (+/-120 nt): {r['nearby_atgs_w_location']}\\n")
        out.write(f"  Upstream(<=300 nt) GC={r['upstream_gc_pct']}%, CpG_OE={r['upstream_cpg_oe']}, TATAAA={r['upstream_tataaa_hits']}\\n\\n")

print(df[['copy','augustus_start','augustus_start_kozak_score_0to2','augustus_sampled_tx_score','best_exon_atg','upstream_gc_pct','upstream_tataaa_hits']].to_string(index=False))
print('Wrote start_site_assessment.tsv and start_site_assessment_summary.txt')
