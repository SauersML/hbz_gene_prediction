#!/usr/bin/env python3
import os, re
from pathlib import Path
from alphagenome.models import dna_client

FA='/Users/user/Downloads/hbz/HBZ_region_UP005hap1.fasta'

def read_fa(p):
    s=[]
    for l in Path(p).read_text().splitlines():
        if l.startswith('>'): continue
        s.append(re.sub(r'[^ACGTNacgtn]','',l))
    return ''.join(s).upper()

seq=read_fa(FA)
target=dna_client.SEQUENCE_LENGTH_100KB
left=(target-len(seq))//2
seq='N'*left + seq + 'N'*(target-len(seq)-left)

terms_list=[
 ['UBERON:0000178'], # blood
 ['UBERON:0002371'], # bone marrow
 ['UBERON:0001114'], # right liver lobe
 ['UBERON:0002106'], # spleen
 ['UBERON:0001987'], # placenta
 ['CL:0000765'],     # erythroblast
 ['CL:0000232'],     # erythrocyte
 ['CL:0000037'],     # hematopoietic stem cell
 ['UBERON:0000955'], # brain
 ['UBERON:0002048'], # lung
]

key=os.getenv('ALPHAGENOME_API_KEY')
if not key:
    raise SystemExit('need ALPHAGENOME_API_KEY')
model=dna_client.create(key)

for terms in terms_list:
    try:
        out=model.predict_sequence(sequence=seq,requested_outputs=[dna_client.OutputType.RNA_SEQ,dna_client.OutputType.CAGE,dna_client.OutputType.PROCAP,dna_client.OutputType.SPLICE_SITES],ontology_terms=terms)
        def shape(obj_name):
            if not hasattr(out,obj_name): return None
            td=getattr(out,obj_name)
            return tuple(td.values.shape)
        print(f"{terms[0]}\trna_seq={shape('rna_seq')}\tcage={shape('cage')}\tprocap={shape('procap')}\tsplice_sites={shape('splice_sites')}")
    except Exception as e:
        print(f"{terms[0]}\tERROR\t{str(e).replace(chr(10),' ')[:220]}")
