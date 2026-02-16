#!/usr/bin/env python3
import re
from itertools import combinations

STOP = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1].upper()


def read_fasta(path):
    out = []
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name:
                    out.append((name, ''.join(seq).upper()))
                name = line[1:].strip()
                seq = []
            else:
                seq.append(re.sub(r'[^ACGTNacgtn]', '', line))
    if name:
        out.append((name, ''.join(seq).upper()))
    return out


def find_orfs(seq, min_aa=30):
    res = []
    n = len(seq)
    for frame in range(3):
        i = frame
        while i <= n - 3:
            codon = seq[i:i+3]
            if codon == 'ATG':
                j = i + 3
                while j <= n - 3:
                    c = seq[j:j+3]
                    if c in STOP:
                        aa = (j + 3 - i) // 3
                        if aa >= min_aa:
                            res.append((i, j + 3, frame, aa))
                        break
                    j += 3
            i += 3
    res.sort(key=lambda x: x[3], reverse=True)
    return res


def canonical_introns(seq, min_len=40, max_len=1200):
    # donor GT, acceptor AG
    n = len(seq)
    donors = [m.start() for m in re.finditer('GT', seq)]
    acceptors = [m.start() for m in re.finditer('AG', seq)]
    out = []
    for d in donors:
        for a in acceptors:
            if a <= d + 1:
                continue
            L = (a + 2) - d
            if min_len <= L <= max_len:
                out.append((d, a + 2))
    return out


def splice_seq(seq, introns):
    introns = sorted(introns)
    pieces = []
    cur = 0
    for s, e in introns:
        pieces.append(seq[cur:s])
        cur = e
    pieces.append(seq[cur:])
    return ''.join(pieces)


def lift_pos(spliced_pos, introns):
    # map 0-based position in spliced seq back to genomic coordinate
    shift = 0
    cur_sp = 0
    for s, e in sorted(introns):
        exon_len = s - (cur_sp + shift)
        if spliced_pos < cur_sp + exon_len:
            return spliced_pos + shift
        cur_sp += exon_len
        shift += (e - s)
    return spliced_pos + shift


def best_spliced_orf(seq, min_aa=30, max_introns=2):
    intr = canonical_introns(seq)
    # heuristic prefilter: introns whose boundaries preserve some frame for coding starts
    # keep only top plausible by size to limit combinatorics
    intr = sorted(intr, key=lambda x: x[1]-x[0])
    if len(intr) > 500:
        intr = intr[:500]

    best = None
    # no intron case
    for k in range(0, max_introns + 1):
        if k == 0:
            combos = [()]
        else:
            combos = combinations(intr, k)
        tested = 0
        for combo in combos:
            combo = sorted(combo)
            ok = True
            for i in range(1, len(combo)):
                if combo[i][0] < combo[i-1][1]:
                    ok = False
                    break
            if not ok:
                continue
            sp = splice_seq(seq, combo)
            orfs = find_orfs(sp, min_aa=min_aa)
            tested += 1
            if not orfs:
                continue
            top = orfs[0]
            score = top[3]
            if (best is None) or (score > best['aa']):
                s, e, fr, aa = top
                best = {
                    'aa': aa,
                    'spliced_start': s,
                    'spliced_end': e,
                    'frame': fr,
                    'introns': combo,
                    'tested': tested,
                    'spliced_len': len(sp),
                    'gen_start': lift_pos(s, combo),
                    'gen_end': lift_pos(e-1, combo)+1,
                }
        # keep runtime contained if huge
        if best and best['aa'] >= 100:
            break
    return best


def summarize_file(path):
    recs = read_fasta(path)
    rows = []
    for name, seq in recs:
        fwd_orfs = find_orfs(seq, min_aa=20)
        rev_orfs = find_orfs(revcomp(seq), min_aa=20)
        best_fwd = fwd_orfs[0] if fwd_orfs else None
        best_rev = rev_orfs[0] if rev_orfs else None

        sp_fwd = best_spliced_orf(seq, min_aa=20, max_introns=2)
        sp_rev = best_spliced_orf(revcomp(seq), min_aa=20, max_introns=2)

        rows.append((name, len(seq), best_fwd, best_rev, sp_fwd, sp_rev))
    return rows


def format_orf(orf):
    if not orf:
        return 'none'
    s,e,fr,aa = orf
    return f"start={s+1}, end={e}, frame={fr}, aa={aa}"


def format_sp(sp, strand):
    if not sp:
        return 'none'
    intr = ', '.join([f"{s+1}-{e}" for s,e in sp['introns']]) if sp['introns'] else 'none'
    return (
        f"strand={strand}, aa={sp['aa']}, genomic={sp['gen_start']+1}-{sp['gen_end']}, "
        f"spliced={sp['spliced_start']+1}-{sp['spliced_end']}, introns={intr}"
    )


def main():
    import glob
    files = sorted(glob.glob('*.fa')) + sorted(glob.glob('*.fasta'))
    print('file\tlen\tbest_unspliced_fwd\tbest_unspliced_rev\tbest_spliced_fwd\tbest_spliced_rev')
    for f in files:
        for name, n, bf, br, sf, sr in summarize_file(f):
            print(
                f"{f}\t{n}\t{format_orf(bf)}\t{format_orf(br)}\t{format_sp(sf,'+')}\t{format_sp(sr,'-')}"
            )

if __name__ == '__main__':
    main()
