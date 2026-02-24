#!/usr/bin/env python3
import re
import sys
from pathlib import Path

COMP = str.maketrans('ACGTNacgtn','TGCANtgcan')
STOP = {'TAA','TAG','TGA'}
CODON_TABLE = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def read_fa(path):
    name=None;seq=[]
    for line in Path(path).read_text().splitlines():
        if not line: continue
        if line.startswith('>'):
            name=line[1:].strip(); continue
        seq.append(re.sub(r'[^ACGTNacgtn]','',line))
    return name,''.join(seq).upper()

def rc(s):
    return s.translate(COMP)[::-1]

def translate(s):
    aa=[]
    for i in range(0,len(s)-2,3):
        aa.append(CODON_TABLE.get(s[i:i+3],'X'))
    return ''.join(aa)

print('file\tstrand\tgene_start\tgene_end\tintrons\tspliced_cds_nt\torf_status\taa_len\tprotein_prefix')
for gff in sorted(Path('.').glob('*.augustus.gff3')):
    fasta = gff.name.replace('.augustus.gff3','.fa')
    if not Path(fasta).exists():
        print(f"Skipping {gff.name}: missing FASTA {fasta}", file=sys.stderr)
        continue
    _, seq = read_fa(fasta)
    strand = '+'
    cds=[]
    gene_start=gene_end=None
    for line in gff.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split('\t')
        feat=f[2]
        s,e = int(f[3]), int(f[4])
        if feat=='gene':
            gene_start, gene_end = s,e
            strand=f[6]
        if feat=='CDS':
            cds.append((s,e))
    cds.sort(key=lambda x:x[0])
    if strand == '-':
        cds = cds[::-1]

    # introns from cds chain on genomic + coords
    cds_sorted = sorted(cds, key=lambda x:x[0])
    intr=[]
    for i in range(len(cds_sorted)-1):
        intr.append((cds_sorted[i][1]+1, cds_sorted[i+1][0]-1))

    parts=[]
    for s,e in cds:
        frag = seq[s-1:e]
        parts.append(frag)
    spliced = ''.join(parts)
    if strand=='-':
        spliced = rc(spliced)
    aa = translate(spliced)
    status=[]
    if spliced.startswith('ATG'): status.append('starts_ATG')
    if spliced[-3:] in STOP: status.append('ends_STOP')
    internal = aa[:-1].count('*') if aa else 0
    if internal==0: status.append('no_internal_stop')
    status=';'.join(status)
    intr_s = 'none' if not intr else ','.join([f'{a}-{b}' for a,b in intr])
    print(f"{gff.name}\t{strand}\t{gene_start}\t{gene_end}\t{intr_s}\t{len(spliced)}\t{status}\t{len(aa)-1 if aa.endswith('*') else len(aa)}\t{aa[:25]}")
