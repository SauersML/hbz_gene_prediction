# hbz_gene_prediction

Comprehensive analysis of three HBZ-region DNA sequences (`copy1`, `copy2`, `copy3`) using multiple gene prediction approaches, splicing-aware ORF checks, and dotplot alignment.

## 1) Scope and objective
This project answers the following biological questions for three FASTA inputs:

1. Are there gene-like coding models in these DNA segments?
2. If introns are spliced, do valid open reading frames (ORFs) remain?
3. Do independent predictors agree on structure?
4. How similar are `copy1` and `copy2` at the nucleotide level (dotplot)?

The analysis was performed interactively in terminal with reproducible scripts and saved outputs.

## 2) Input data
Primary input FASTA files:

- `copy1.fa` (1930 nt)
- `copy2.fa` (2176 nt)
- `copy3.fa` (2182 nt)

Headers indicate all are from the same broader HBZ region context and are plus-strand slices.

## 3) Analysis chronology (what was done, end-to-end)
This section summarizes the full workflow from the session, in order.

### Phase A: Initial ORF and splice-aware scan
A custom script (`hbz_orf_splice_scan.py`) was created to:

- scan 6 frames for unspliced ORFs,
- enumerate canonical `GT...AG` introns,
- test spliced sequence combinations (up to 2 introns),
- report best ORF in forward/reverse orientations.

Output:

- `hbz_orf_splice_report.tsv`

Initial signal:

- strong long ORFs were detected in all copies,
- splicing did not eliminate coding potential.

### Phase B: Install and run eukaryotic gene predictor
`AUGUSTUS` was installed and run with human parameters:

- command style: `augustus --species=human --gff3=on copyX.fa`

Outputs:

- `copy1.augustus.gff3`
- `copy2.augustus.gff3`
- `copy3.augustus.gff3`

Findings:

- 1 gene/transcript per copy,
- 3 CDS exons + 2 introns per copy,
- coherent start/stop and coding structure in all three.

### Phase C: Spliced CDS extraction and ORF validation
A script (`extract_spliced_orf_from_augustus.py`) extracted AUGUSTUS CDS chains and verified ORF properties.

Output:

- `hbz_spliced_orf_summary.tsv`

Validation checks (all copies):

- starts with `ATG`,
- ends with stop codon,
- no internal stop codons.

### Phase D: Additional predictors (stronger multi-model comparison)
To avoid single-tool bias, two additional predictors were run:

- **GlimmerHMM** (human trained directory)
- **Prodigal** (`-p meta`, prokaryotic-biased reference call set)

Outputs:

- `copy1.glimmer.gff3`, `copy2.glimmer.gff3`, `copy3.glimmer.gff3`
- `copy1.prodigal.gff`, `copy2.prodigal.gff`, `copy3.prodigal.gff`
- peptide/CDS FASTAs from Prodigal

Consensus behavior:

- AUGUSTUS and GlimmerHMM gave near-identical multi-exon models for each copy.
- Prodigal produced mostly single-exon fragmented calls (expected because it is optimized for prokaryotic gene structure, not intron-rich eukaryotic models).

### Phase E: Figure generation
Scripts created publication-style summary figures:

- `make_hbz_figures.py`
- `make_multi_predictor_gene_model_figure.py`

Generated figures:

- `hbz_gene_models.png` (AUGUSTUS model view)
- `hbz_orf_lengths.png` (predicted protein lengths)
- `hbz_all_predictors_gene_models.png` (single consolidated predictor-by-copy panel)
- corresponding `.svg` where available

### Phase F: AUGUSTUS protein FASTA export
AUGUSTUS protein sequences were exported into:

- `augustus_predicted_proteins.fa`

Contains:

- `copy1` predicted protein: 142 aa
- `copy2` predicted protein: 129 aa
- `copy3` predicted protein: 129 aa

### Phase G: Dotplot alignment (`copy1` vs `copy2`)
A k-mer dotplot script was created:

- `dotplot_copy1_vs_copy2.py`

Outputs:

- `copy1_vs_copy2_dotplot.png/.svg`
- `copy1_vs_copy2_dotplot_augustus_exons.png/.svg` (with exon overlays)
- `copy1_vs_copy2_dotplot_stats.txt`

Observed pattern:

- strong forward orientation similarity,
- no reverse-complement signal (no inversion-like relationship),
- multi-diagonal structure consistent with shared blocks plus repeat/offset-rich segments.

## 4) Core quantitative results

### 4.1 AUGUSTUS/GlimmerHMM spliced models

- `copy1`: gene `267-1827`, CDS exons at `267-361`, `1250-1454`, `1699-1827`
  - introns: `362-1249`, `1455-1698`
  - spliced CDS: 429 nt
  - protein: 142 aa

- `copy2`: gene `228-2073`, CDS exons at `228-283`, `1550-1754`, `1945-2073`
  - introns: `284-1549`, `1755-1944`
  - spliced CDS: 390 nt
  - protein: 129 aa

- `copy3`: gene `228-2139`, CDS exons at `228-283`, `1536-1740`, `2011-2139`
  - introns: `284-1535`, `1741-2010`
  - spliced CDS: 390 nt
  - protein: 129 aa

### 4.2 Protein comparison

- `copy2` and `copy3` are identical in predicted AUGUSTUS protein sequence.
- `copy1` shares the same core sequence and is longer at the N-terminus.

## 5) Biological interpretation

1. All three sequences contain robust **protein-coding candidates** when evaluated with eukaryotic predictors and splice-aware ORF validation.
2. Independent eukaryotic models (AUGUSTUS + GlimmerHMM) converge on the same architecture per copy, increasing confidence.
3. `copy2` and `copy3` are highly similar coding forms; `copy1` likely represents an extended/alternate 5' coding configuration.
4. Dotplot evidence supports same-orientation homology between `copy1` and `copy2`.

## 6) Caveats and limitations

- Predictions are **ab initio** (no RNA-seq/cDNA/protein-alignment hints supplied during prediction).
- Results indicate strong gene candidacy, but not final reference-grade annotation.
- Prodigal disagreement is expected due to model domain mismatch (prokaryotic assumptions).

## 7) Reproducibility quickstart
Run from repository root:

```bash
# Splice-aware ORF scan
python3 hbz_orf_splice_scan.py > hbz_orf_splice_report.tsv

# AUGUSTUS models
augustus --species=human --gff3=on copy1.fa > copy1.augustus.gff3
augustus --species=human --gff3=on copy2.fa > copy2.augustus.gff3
augustus --species=human --gff3=on copy3.fa > copy3.augustus.gff3

# GlimmerHMM (path may vary by install)
glimmerhmm copy1.fa /opt/homebrew/Caskroom/miniforge/base/share/glimmerhmm/trained_dir/human -g -o copy1.glimmer.gff3
glimmerhmm copy2.fa /opt/homebrew/Caskroom/miniforge/base/share/glimmerhmm/trained_dir/human -g -o copy2.glimmer.gff3
glimmerhmm copy3.fa /opt/homebrew/Caskroom/miniforge/base/share/glimmerhmm/trained_dir/human -g -o copy3.glimmer.gff3

# Prodigal
prodigal -i copy1.fa -f gff -o copy1.prodigal.gff -a copy1.prodigal.aa.fa -d copy1.prodigal.cds.fa -p meta
prodigal -i copy2.fa -f gff -o copy2.prodigal.gff -a copy2.prodigal.aa.fa -d copy2.prodigal.cds.fa -p meta
prodigal -i copy3.fa -f gff -o copy3.prodigal.gff -a copy3.prodigal.aa.fa -d copy3.prodigal.cds.fa -p meta

# Build figures/tables
python3 make_hbz_figures.py
python3 make_multi_predictor_gene_model_figure.py
python3 dotplot_copy1_vs_copy2.py
```

## 8) File guide

### Inputs
- `copy1.fa`, `copy2.fa`, `copy3.fa`

### Predictor outputs
- `copy*.augustus.gff3`
- `copy*.glimmer.gff3`
- `copy*.prodigal.gff`
- `copy*.prodigal.aa.fa`
- `copy*.prodigal.cds.fa`

### Summaries/tables
- `hbz_gene_prediction_summary.tsv`
- `hbz_multi_predictor_summary.tsv`
- `hbz_spliced_orf_summary.tsv`
- `hbz_orf_splice_report.tsv`
- `copy1_vs_copy2_dotplot_stats.txt`

### Figures
- `hbz_gene_models.png`
- `hbz_orf_lengths.png`
- `hbz_all_predictors_gene_models.png`
- `copy1_vs_copy2_dotplot.png`
- `copy1_vs_copy2_dotplot_augustus_exons.png`

### Scripts
- `hbz_orf_splice_scan.py`
- `extract_spliced_orf_from_augustus.py`
- `compare_augustus_proteins.py`
- `make_hbz_figures.py`
- `make_multi_predictor_gene_model_figure.py`
- `dotplot_copy1_vs_copy2.py`

## 9) Bottom line
All three HBZ-region copies show strong coding potential with consistent eukaryotic multi-exon gene models. The most stable conclusion from this analysis is a shared spliced coding architecture across copies, with `copy1` encoding a longer predicted protein variant than `copy2/3`.

## 10) Start-codon / initiation plausibility follow-up
To address whether AUGUSTUS starts are biologically plausible for initiation, an exon-aware start assessment was added.

New files:

- `start_site_assessment.py`
- `start_site_assessment.tsv`
- `start_site_assessment_summary.txt`
- `start_site_assessment.png`
- `start_site_assessment.svg`
- `copy*.augustus.sampled.gff3` (AUGUSTUS sampling runs)

What this follow-up checks:

1. AUGUSTUS sampled transcript score (`--alternatives-from-sampling=true --sample=200`)
2. Kozak core score at AUGUSTUS start (`-3 in {A,G}` and `+4 == G`)
3. Alternative ATGs restricted to **predicted CDS exons** (mRNA-feasible starts)
4. Nearby ATGs annotated as exon/intron/other to avoid intron false leads
5. Upstream promoter proxies (GC%, CpG observed/expected, TATAAA count)

Key result:

- `copy2/3` have apparently stronger nearby ATGs in genomic DNA (e.g., around position 331), but those fall in **predicted introns**, so they are not start codons in the mature spliced transcript.
- Exon-contained ATG options in `copy2/3` collapse to the AUGUSTUS start at 228.
- AUGUSTUS still returns the same single model in sampling mode (`~0.89-0.905` score range), indicating internal model stability despite weak Kozak context.

Interpretation:

- Current sequence-only evidence supports a stable **gene model** and a feasible coding start under that model.
- Sequence-only evidence does **not** prove real transcription start site usage; expression/TSS assays are required for that.
