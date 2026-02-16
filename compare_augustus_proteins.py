#!/usr/bin/env python3
import re
from pathlib import Path

# get protein sequence from AUGUSTUS comment block
prots={}
for gff in sorted(Path('.').glob('copy*.augustus.gff3')):
    seq=[]
    cap=False
    for line in gff.read_text().splitlines():
        if line.startswith('# protein sequence = ['):
            cap=True
            part=line.split('[',1)[1]
            if ']' in part:
                seq.append(part.split(']')[0].replace(' ',''))
                cap=False
            else:
                seq.append(part.replace(' ',''))
            continue
        if cap:
            text=line[1:].strip() if line.startswith('#') else line.strip()
            if ']' in text:
                seq.append(text.split(']')[0].replace(' ',''))
                cap=False
            else:
                seq.append(text.replace(' ',''))
    prots[gff.stem.replace('.augustus','')]=''.join(seq)

names=sorted(prots)
print('sample\taa_len')
for n in names:
    print(f'{n}\t{len(prots[n])}')

print('\nPairwise identity over overlap')
for i,a in enumerate(names):
    for b in names[i+1:]:
        s1,s2=prots[a],prots[b]
        m=min(len(s1),len(s2))
        ident=sum(1 for x,y in zip(s1[:m],s2[:m]) if x==y)
        pid=100*ident/m if m else 0
        print(f'{a} vs {b}: {pid:.1f}% over {m} aa (len {len(s1)} vs {len(s2)})')
