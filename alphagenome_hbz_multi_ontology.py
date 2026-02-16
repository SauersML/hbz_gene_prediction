#!/usr/bin/env python3
import os, re, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from alphagenome.models import dna_client

FA='/Users/user/Downloads/hbz/HBZ_region_UP005hap1.fasta'
OUT=Path('alphagenome_hbz_multi_ontology')
OUT.mkdir(exist_ok=True)

HBZ_WINDOWS=[('copy1_region',1705,3635),('copy2_region',12466,14642),('copy3_region',22484,24666)]
ONTOLOGIES=[
    'UBERON:0001987', # placenta
    'UBERON:0002371', # bone marrow
    'UBERON:0002106', # spleen
    'UBERON:0001114', # right liver lobe
    'UBERON:0000178', # blood
    'UBERON:0002048', # lung
    'UBERON:0000955', # brain
]
OUTPUTS=[dna_client.OutputType.RNA_SEQ,dna_client.OutputType.CAGE,dna_client.OutputType.SPLICE_SITES]


def read_fa(p):
    s=[]
    for l in Path(p).read_text().splitlines():
        if l.startswith('>'): continue
        s.append(re.sub(r'[^ACGTNacgtn]','',l))
    return ''.join(s).upper()

def pad(seq, target):
    left=(target-len(seq))//2
    return 'N'*left+seq+'N'*(target-len(seq)-left), left

seq=read_fa(FA)
orig_len=len(seq)
model_len=dna_client.SEQUENCE_LENGTH_100KB if orig_len<=dna_client.SEQUENCE_LENGTH_100KB else dna_client.SEQUENCE_LENGTH_500KB
seq_pad,pad_left=pad(seq,model_len)

key=os.getenv('ALPHAGENOME_API_KEY')
if not key:
    raise SystemExit('Missing ALPHAGENOME_API_KEY')
model=dna_client.create(key)

rows=[]
fail=[]
for ont in ONTOLOGIES:
    try:
        out=model.predict_sequence(sequence=seq_pad, requested_outputs=OUTPUTS, ontology_terms=[ont])
    except Exception as e:
        fail.append({'ontology':ont,'error':str(e)})
        continue

    for tname in ['rna_seq','cage','splice_sites']:
        if not hasattr(out,tname):
            continue
        td=getattr(out,tname)
        vals=td.values[pad_left:pad_left+orig_len,:]
        if vals.shape[1]==0:
            for w,s,e in HBZ_WINDOWS:
                rows.append({'ontology':ont,'track':tname,'window':w,'start':s,'end':e,'window_mean':np.nan,'global_mean':np.nan,'enrichment':np.nan,'tracks_n':0})
            continue

        per=vals.mean(axis=1)
        gmean=float(np.mean(per)) if per.size else np.nan
        for w,s,e in HBZ_WINDOWS:
            ww=per[s-1:e]
            wmean=float(np.mean(ww)) if ww.size else np.nan
            enrich=(wmean/gmean) if (gmean and not np.isnan(gmean)) else np.nan
            rows.append({'ontology':ont,'track':tname,'window':w,'start':s,'end':e,'window_mean':wmean,'global_mean':gmean,'enrichment':enrich,'tracks_n':vals.shape[1]})

pd.DataFrame(rows).to_csv(OUT/'window_enrichment_by_ontology.tsv',sep='\t',index=False)
Path(OUT/'failed_ontologies.json').write_text(json.dumps(fail,indent=2))

# heuristic calls
calls=[]
df=pd.DataFrame(rows)
for ont in df['ontology'].dropna().unique():
    for w,_,_ in HBZ_WINDOWS:
        sub=df[(df.ontology==ont)&(df.window==w)]
        rna=float(sub[sub.track=='rna_seq']['enrichment'].iloc[0]) if not sub[sub.track=='rna_seq'].empty else np.nan
        cage=float(sub[sub.track=='cage']['enrichment'].iloc[0]) if not sub[sub.track=='cage'].empty else np.nan
        splice=float(sub[sub.track=='splice_sites']['enrichment'].iloc[0]) if not sub[sub.track=='splice_sites'].empty else np.nan
        flags={
            'rna_high': (not np.isnan(rna)) and rna>=1.5,
            'cage_high': (not np.isnan(cage)) and cage>=1.5,
            'splice_high': (not np.isnan(splice)) and splice>=1.5,
        }
        pos=sum(flags.values())
        verdict='supportive' if pos>=2 else ('weak' if pos==1 else 'not_supportive')
        calls.append({'ontology':ont,'window':w,'rna_enrich':rna,'cage_enrich':cage,'splice_enrich':splice,**flags,'verdict':verdict})

cdf=pd.DataFrame(calls)
cdf.to_csv(OUT/'transcription_support_calls.tsv',sep='\t',index=False)

# heatmap figure for RNA and CAGE
plot_df=df[df['track'].isin(['rna_seq','cage'])].copy()
plot_df['key']=plot_df['ontology']+' | '+plot_df['track']
mat=plot_df.pivot_table(index='key',columns='window',values='enrichment',aggfunc='mean')
plt.figure(figsize=(9,6),dpi=220)
sns.heatmap(mat,annot=True,fmt='.2f',cmap='RdYlGn',center=1.0,vmin=0,vmax=max(2.5,np.nanmax(mat.values)))
plt.title('AlphaGenome HBZ window enrichment by ontology (RNA/CAGE)')
plt.tight_layout()
plt.savefig(OUT/'hbz_transcription_enrichment_heatmap.png')
plt.savefig(OUT/'hbz_transcription_enrichment_heatmap.svg')

print('done')
print('rows',len(rows),'calls',len(calls),'fails',len(fail))
print((OUT/'transcription_support_calls.tsv').as_posix())
